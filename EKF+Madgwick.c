#include <assert.h>
#include "KalmanMadgwick.h"
#include "Matrix.h"
#include "Kalman.h"
#include "params_noise.h"
#include "MadgwickAHRS.h"

KalmanMadgwickFilter_t *KalmanMadgwickAlloc(double x,
                                            double y,
                                            double xVel,
                                            double yVel,
                                            double accDev,
                                            double posDev,
                                            double timeStamp)
{
    KalmanMadgwickFilter_t *f = (KalmanMadgwickFilter_t *)malloc(sizeof(KalmanMadgwickFilter_t));
    assert(f);
    f->Rot = MatrixAlloc(3, 3); // rotation matrix
    f->E = MatrixAlloc(3, 3);   // rotation matrix
    f->kf = KalmanFilterCreate(15, 3, 0);
    /*initialization*/
    f->predictTime = f->updateTime = timeStamp;
    f->predictCount = 0;
    f->accDev = accDev;

    MatrixSet(f->kf->Xk_k,
              x, y, xVel, yVel);
    MatrixSetIdentityDiag(f->kf->H); // state has 4d and measurement has 4d too. so here is identity

    MatrixSetIdentity(f->kf->Pk_k);

    f->kf->Pk_k->data[0][0] = 0.5;
    f->kf->Pk_k->data[1][1] = 0.5;
    f->kf->Pk_k->data[2][2] = 0.5;

    f->kf->Pk_k->data[6][6] = 10;
    f->kf->Pk_k->data[7][7] = 10;
    f->kf->Pk_k->data[8][8] = 10;

    f->kf->Pk_k->data[9 ][9 ] = 500;
    f->kf->Pk_k->data[10][10] = 500;
    f->kf->Pk_k->data[11][11] = 500;

    f->kf->Pk_k->data[12][12] = 500;
    f->kf->Pk_k->data[13][13] = 500;
    f->kf->Pk_k->data[14][14] = 500;

    return f;
}
//////////////////////////////////////////////////////////////////////////

void KalmanMadgwickFree(KalmanMadgwickFilter_t *k)
{
    assert(k);
    KalmanFilterFree(k->kf);
    free(k);
}

//////////////////////////////////////////////////////////////////////////

static void rebuildRot(KalmanMadgwickFilter_t *f)
{
    MatrixSet(f->Rot,
              2 * pow(q0, 2) - 1 + 2 * pow(q1, 2), 2 * (q1 * q2 + q0 * q3), 2 * (q1 * q3 - q0 * q2),
              2 * (q1 * q2 - q0 * q3), 2 * pow(q0, 2) - 1 + 2 * pow(q2, 2), 2 * (q2 * q3 + q0 * q1),
              2 * (q1 * q3 + q0 * q2), 2 * (q2 * q3 - q0 * q1), 2 * pow(q0, 2) - 1 + 2 * pow(q3, 2));
}
//////////////////////////////////////////////////////////////////////////

static void rebuildE(KalmanMadgwickFilter_t *f)
{
    double phi0 = atan2(f->Rot->data[2][1], f->Rot->data[2][2]);
    double phi1 = -atan(f->Rot->data[2][0] / sqrt(1 - pow(f->Rot->data[2][0], 2)));
    double phi2 = atan2(f->Rot->data[1][0], f->Rot->data[0][0]);

    MatrixSet(f->E,
              1.0, sin(phi0) * tan(phi1), cos(phi0) * tan(phi1),
              0.0, cos(phi0), -1 * sin(phi0),
              0.0, sin(phi0) / cos(phi1), cos(phi0) / cos(phi1));
}
//////////////////////////////////////////////////////////////////////////

static void rebuildF(KalmanMadgwickFilter_t *f, double dt)
{
    rebuildRot(f);
    rebuildE(f);

    MatrixSetIdentityDiag(f->kf->F);

    f->kf->F->data[0][3] = dt;
    f->kf->F->data[1][4] = dt;
    f->kf->F->data[2][5] = dt;

    MatrixScale(f->Rot, (-1 * dt));

    f->kf->F->data[3][12] = f->Rot->data[0][0];
    f->kf->F->data[3][13] = f->Rot->data[0][1];
    f->kf->F->data[3][14] = f->Rot->data[0][2];

    f->kf->F->data[4][12] = f->Rot->data[1][0];
    f->kf->F->data[4][13] = f->Rot->data[1][1];
    f->kf->F->data[4][14] = f->Rot->data[1][2];

    f->kf->F->data[5][12] = f->Rot->data[2][0];
    f->kf->F->data[5][13] = f->Rot->data[2][1];
    f->kf->F->data[5][14] = f->Rot->data[2][2];

    rebuildRot(f);
    MatrixScale(f->Rot, (-1 * dt * dt));
    MatrixScale(f->Rot, 0.5);

    f->kf->F->data[0][12] = f->Rot->data[0][0];
    f->kf->F->data[0][13] = f->Rot->data[0][1];
    f->kf->F->data[0][14] = f->Rot->data[0][2];

    f->kf->F->data[1][12] = f->Rot->data[1][0];
    f->kf->F->data[1][13] = f->Rot->data[1][1];
    f->kf->F->data[1][14] = f->Rot->data[1][2];

    f->kf->F->data[2][12] = f->Rot->data[2][0];
    f->kf->F->data[2][13] = f->Rot->data[2][1];
    f->kf->F->data[2][14] = f->Rot->data[2][2];

    MatrixScale(f->E, (-1 * dt));
    f->kf->F->data[6][9] = f->E->data[0][0];
    f->kf->F->data[6][10] = f->E->data[0][1];
    f->kf->F->data[6][11] = f->E->data[0][2];

    f->kf->F->data[7][9] = f->E->data[1][0];
    f->kf->F->data[7][10] = f->E->data[1][1];
    f->kf->F->data[7][11] = f->E->data[1][2];

    f->kf->F->data[8][9] = f->E->data[2][0];
    f->kf->F->data[8][10] = f->E->data[2][1];
    f->kf->F->data[8][11] = f->E->data[2][2];
}
//////////////////////////////////////////////////////////////////////////

static void rebuildU(KalmanMadgwickFilter_t *f,
                     double xAcc,
                     double yAcc)
{
    MatrixSet(f->kf->Uk,
              xAcc,
              yAcc);
}
//////////////////////////////////////////////////////////////////////////

static void rebuildB(KalmanMadgwickFilter_t *f,
                     double dt)
{
    double dt2 = 0.5 * dt * dt;
    MatrixSet(f->kf->B,
              dt2, 0.0,
              0.0, dt2,
              dt, 0.0,
              0.0, dt);
}
//////////////////////////////////////////////////////////////////////////

static void rebuildR(KalmanMadgwickFilter_t *f,
                     double posSigma)
{
    //  MatrixSetIdentity(f->kf->R);
    //  MatrixScale(f->kf->R, posSigma*posSigma);
    double velSigma = posSigma * 1.0e-01;
    MatrixSet(f->kf->R,
              pow(RTK_NOISE[X], 2), 0.0, 0.0,
              0.0, pow(RTK_NOISE[Y], 2), 0.0,
              0.0, 0.0, pow(RTK_NOISE[Z], 2));
}
//////////////////////////////////////////////////////////////////////////

static void rebuildQ(KalmanMadgwickFilter_t *f)
{
    f->kf->Q->data[3][3] = g2ms2(NOISE_ACC[X]);
    f->kf->Q->data[4][4] = g2ms2(NOISE_ACC[Y]);
    f->kf->Q->data[5][5] = g2ms2(NOISE_ACC[Z]);

    f->kf->Q->data[6][6] = Degree2Rad(NOISE_GYR[X]);
    f->kf->Q->data[7][7] = Degree2Rad(NOISE_GYR[Y]);
    f->kf->Q->data[8][8] = Degree2Rad(NOISE_GYR[Z]);

    f->kf->Q->data[9][9] = Degree2Rad(BIAS_GYR[X]);
    f->kf->Q->data[10][10] = Degree2Rad(BIAS_GYR[Y]);
    f->kf->Q->data[11][11] = Degree2Rad(BIAS_GYR[Z]);

    f->kf->Q->data[12][12] = g2ms2(BIAS_ACC[X]);
    f->kf->Q->data[13][13] = g2ms2(BIAS_ACC[Y]);
    f->kf->Q->data[14][14] = g2ms2(BIAS_ACC[Z]);
}

void KalmanMadgwickPredict(KalmanMadgwickFilter_t *k,
                           double timeNow,
                           double xAcc,
                           double yAcc)
{
    double dt = (timeNow - k->predictTime) / 1.0e+3;

    rebuildF(k, dt);
    rebuildB(k, dt);
    rebuildU(k, xAcc, yAcc);

    ++k->predictCount;
    rebuildQ(k);

    k->predictTime = timeNow;
    KalmanFilterPredict(k->kf);
    MatrixCopy(k->kf->Xk_km1, k->kf->Xk_k);
}
//////////////////////////////////////////////////////////////////////////

void KalmanMadgwickUpdate(KalmanMadgwickFilter_t *k,
                          double timeNow,
                          double x,
                          double y,
                          double xVel,
                          double yVel,
                          double posDev)
{
    //  double dt = timeNow - k->updateTime;
    k->predictCount = 0;
    k->updateTime = timeNow;
    rebuildR(k, posDev);
    MatrixSet(k->kf->Zk, x, y, xVel, yVel);
    KalmanFilterUpdate(k->kf);
}
//////////////////////////////////////////////////////////////////////////

double KalmanMadgwickGetX(const KalmanMadgwickFilter_t *k)
{
    return k->kf->Xk_k->data[0][0];
}

double KalmanMadgwickGetY(const KalmanMadgwickFilter_t *k)
{
    return k->kf->Xk_k->data[1][0];
}

double KalmanMadgwickGetXVel(const KalmanMadgwickFilter_t *k)
{
    return k->kf->Xk_k->data[2][0];
}

double KalmanMadgwickGetYVel(const KalmanMadgwickFilter_t *k)
{
    return k->kf->Xk_k->data[3][0];
}
//////////////////////////////////////////////////////////////////////////