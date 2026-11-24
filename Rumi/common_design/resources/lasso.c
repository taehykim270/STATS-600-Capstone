#include <R.h>
#include <math.h>
#include <R_ext/Lapack.h>

/* Coordinate descent on active columns, with full KKT checks. G = X'X/n. */
void fit_path(double *g, double *xy, int *size, double *lambda, int *count,
              double *path, double *errors) {
    int p = *size, nactive = 0;
    int *active = (int *) R_alloc(p, sizeof(int));
    int *included = (int *) R_alloc(p, sizeof(int));
    double *b = (double *) R_alloc(p, sizeof(double));
    double *res = (double *) R_alloc(p, sizeof(double));
    /* A restricted Newton step solves the smooth quadratic on the current
       nonzero support. Stop at the first sign change; full KKT checks follow. */
    double *subgram = (double *) R_alloc((size_t)p*p, sizeof(double));
    double *proposal = (double *) R_alloc(p,sizeof(double));
    int *support = (int *) R_alloc(p,sizeof(int));
    for (int j=0; j<p; ++j) { b[j]=0; included[j]=0; res[j]=xy[j]; }
    for (int k=0; k<*count; ++k) {
        double lam=lambda[k], err=0;
        for (int outer=0; outer<2000; ++outer) {
            /* Recompute correlations and bring in at most eight new columns. */
            for (int i=0; i<p; ++i) res[i]=xy[i];
            for (int t=0; t<nactive; ++t) {
                int j=active[t];
                if (b[j]!=0) for (int i=0; i<p; ++i) res[i]-=g[i+p*j]*b[j];
            }
            err=0;
            for (int j=0; j<p; ++j) {
                double v=b[j]==0 ? fmax(fabs(res[j])-lam,0) : fabs(res[j]-copysign(lam,b[j]));
                err=fmax(err,v);
            }
            if (err<=1e-5) break;
            for (int add=0; add<8; ++add) {
                int best=-1; double excess=1e-5;
                for (int j=0; j<p; ++j) if (!included[j] && fabs(res[j])-lam>excess) {
                    best=j; excess=fabs(res[j])-lam;
                }
                if (best<0) break;
                active[nactive++]=best; included[best]=1;
            }
            for (int sweep=0; sweep<1000; ++sweep) {
                double change=0;
                for (int t=0; t<nactive; ++t) {
                    int j=active[t];
                    double z=res[j]+g[j+p*j]*b[j];
                    double next=copysign(fmax(fabs(z)-lam,0),z)/g[j+p*j];
                    double delta=next-b[j];
                    if (delta==0) continue;
                    b[j]=next;
                    change=fmax(change,fabs(delta));
                    for (int u=0; u<nactive; ++u) {
                        int i=active[u]; res[i]-=g[i+p*j]*delta;
                    }
                }
                if (change<1e-8) break;
                if (sweep % 12 == 11) {
                    int n=0, one=1, info=0;
                    for(int t=0;t<nactive;t++) if(b[active[t]]!=0) support[n++]=active[t];
                    if(!n) break;
                    for(int j=0;j<n;j++) {
                        int col=support[j];
                        proposal[j]=xy[col]-copysign(lam,b[col]);
                        for(int i=0;i<=j;i++) subgram[i+n*j]=g[support[i]+p*col];
                    }
                    F77_CALL(dposv)("U",&n,&one,subgram,&n,proposal,&n,&info FCONE);
                    if(info==0) {
                        double step=1;
                        for(int i=0;i<n;i++) {
                            double old=b[support[i]];
                            if(proposal[i]*old<0) step=fmin(step,old/(old-proposal[i]));
                        }
                        for(int i=0;i<n;i++) b[support[i]]+=step*(proposal[i]-b[support[i]]);
                        break;
                    }
                }
            }
            R_CheckUserInterrupt();
        }
        errors[k]=err;
        for (int j=0; j<p; ++j) path[j+p*k]=b[j];
    }
}
