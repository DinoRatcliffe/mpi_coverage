//-----------------------------------------------------------------------------
// c c o v e r a g e   --   calculate the coverage of feature detectors
//-----------------------------------------------------------------------------
// COMPILE: gcc -Wall -g -o ccoverage ccoverage.c -lm
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define MAXFN       200
#define MAXLINE     2000
#define MAXLOCS     200000

#define JOBTAG      0
#define RESULTTAG   1
#define QUITTAG     2
#define ROOTRANK    0

//globals
char *opnames[] = {"ebr", "ibr", "mser", "sfop"};
int nopnames = 4;
char *imsets[] = {"bark", "bikes", "boat", "graf", "leuv", "trees", "ubc",
    "wall"};
int nimsets = 8;

//mpi globals
int my_id, ierr, num_procs, running;
    
int load_locations (char *fn, float locs[][2], int first)
{
    char line[MAXLINE], *kwd;
    int nlocs;
    FILE *f;

    if ((f = fopen (fn, "r")) == NULL) {
        fprintf (stderr, "Cannot open %s!\n", fn);
        exit (EXIT_FAILURE);
    }
    fgets (line, MAXLINE, f);
    fgets (line, MAXLINE, f);
    nlocs = first;
    while (fgets (line, MAXLINE, f) != NULL) {
        kwd = strtok (line, " ");
        locs[nlocs][1] = atof (kwd);
        kwd = strtok ((char *) NULL, " ");
        locs[nlocs][0] = atof (kwd);
        nlocs += 1;
        if (nlocs >= MAXLOCS) {
            fprintf (stderr, "Not enough space to hold all locations: %d!\n", nlocs);
            exit (EXIT_FAILURE);
        }
    }
    fclose (f);
    return nlocs;
}

void coverage (float locs[][2], int nlocs, float *result, int *npaths, int *ncoin)
{
    int np_overall = 0, nc = 0, np, i, j;
    float hmsum = 0.0, dsum, hm, d, y1, x1, y2, x2;

    for (i = 0; i < nlocs; i++) {
        y1 = locs[i][0];
        x1 = locs[i][1];
        np = 0;
        dsum = 0.0;
        for (j = 0; j < nlocs; j++) {
            if (i != j) {
                y2 = locs[j][0] - y1;
                x2 = locs[j][1] - x1;
                d = sqrt (y2 * y2 + x2 * x2);
                if (d <= 0.0) {
                    nc += 1;
                } else {
                    dsum += 1.0 / d;
                    np += 1;
                    np_overall += 1;
                }
            }
        }
        hm = np / dsum;
        hmsum += 1.0 / hm;
    }
    *result = nlocs / hmsum;
    *npaths = np_overall;
    *ncoin = nc;
}

float mean_coverage(int *opindex, int nopindex)
{
    int nfiles = 0, i, j, fno, nlocs, np, nc;
    float locs[MAXLOCS][2];
    float csum = 0.0, cov;
    char fns[nopindex][MAXFN];

    for (i = 0; i < nimsets; i++) {
        for (fno = 1; fno < 7; fno++) {
            int params[7] = {opindex[0], opindex[1], opindex[2], opindex[3], nopindex, i, fno};
            MPI_Request request;
            int proc = (rand() % (num_procs-1)) + 1;
            //printf("sending: %i\n", proc);
            MPI_Isend(&params, 7, MPI_INT, proc, JOBTAG, MPI_COMM_WORLD, &request);
            nfiles++;
        }
    }

    int message_num;
    for(message_num = 0; message_num < nimsets*6; message_num++) {
        MPI_Status status;
        float cov;
        MPI_Recv(&cov, 1, MPI_FLOAT, MPI_ANY_SOURCE, RESULTTAG, MPI_COMM_WORLD, &status);
        //printf("recieved: %i\n", status.MPI_SOURCE);
        csum += cov; 
    }

    return csum / nfiles;
}

void handle_coverage_job(int *opindex, int nopindex, int i, int fno) {
    int nlocs = 0, np, nc, j;
    float cov;
    float locs[MAXLOCS][2];
    char fns[nopindex][MAXFN];

    for (j = 0; j < nopindex; j++) {
       sprintf(fns[j], "data/%s_%s%d.txt", opnames[opindex[j]], imsets[i], fno);
       nlocs = load_locations(fns[j], locs, nlocs);
    }
    coverage(locs, nlocs, &cov, &np, &nc);
            
    MPI_Request request;
    MPI_Isend(&cov, 1, MPI_FLOAT, ROOTRANK, RESULTTAG, MPI_COMM_WORLD, &request);
}

void handle_message() {
    MPI_Status status;
    int params[7];
    MPI_Recv(&params, 7, MPI_INT, ROOTRANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    switch(status.MPI_TAG) {
        case JOBTAG:
            handle_coverage_job(params, params[4], params[5], params[6]);
            break;
        case QUITTAG:
            running = 0;
            break;
    }
}

int main (int argc, char **argv)
{
    int i1, i2, i3, i4;
    float mc;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (my_id == ROOTRANK) {
        // Calculate the coverage for each individual operator.
        for (i1 = 0; i1 < nopnames; i1++) {
            int currentops[1] = {i1};
            mc = mean_coverage (currentops, 1);
            printf ("%s: %f\n", opnames[i1], mc);
        }

        // Calculate the coverage for each combination of pairs of operators.
        for (i1 = 0; i1 < nopnames; i1++) {
            for (i2 = i1+1; i2 < nopnames; i2++) {
                int currentops[2] = {i1, i2};
                mc = mean_coverage (currentops, 2);
                printf ("%s + %s: %f\n", opnames[i1], opnames[i2], mc);
            }
        }

        // Calculate the coverage for each combination of triplets of operators.
        for (i1 = 0; i1 < nopnames; i1++) {
            for (i2 = i1+1; i2 < nopnames; i2++) {
                for (i3 = i2+1; i3 < nopnames; i3++) {
                    int currentops[3] = {i1, i2, i3};
                    mc = mean_coverage (currentops, 3);
                    printf ("%s + %s + %s: %f\n", opnames[i1], opnames[i2], opnames[i3], mc);
                }
            }
        }

        // Calculate the coverage for each combination of quadruplets of operators.
        for (i1 = 0; i1 < nopnames; i1++) {
            for (i2 = i1+1; i2 < nopnames; i2++) {
                for (i3 = i2+1; i3 < nopnames; i3++) {
                    for (i4 = i3+1; i4 < nopnames; i4++) {
                        int currentops[4] = {i1, i2, i3, i4};
                        mc = mean_coverage (currentops, 4);
                        printf ("%s + %s + %s + %s: %f\n", opnames[i1], opnames[i2],
                                opnames[i3], opnames[i4], mc);
                    }
                }
            }
        }

        int current_proc;
        MPI_Request request;
        for (current_proc = 1; current_proc < num_procs; current_proc++) {
            MPI_Isend(&current_proc, 1, MPI_INT, current_proc, QUITTAG, MPI_COMM_WORLD, &request);
        }
    } else {
        running = 1;
        while ( running ) {
            handle_message();
        }
    }

    ierr = MPI_Finalize();

    return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------------
// End of ccoverage.c
//-----------------------------------------------------------------------------
