#include <assert.h>
#include "KalmanMadgwick.h"
#include "Matrix.h"
#include "Kalman.h"
#include "params_noise.h"
#include "MadgwickAHRS.h"

#define STATE_COUNT 15

volatile float phi0 = 0.0f, phi1 = 0.0f, phi2 = 0.0f;

KalmanMadgwickFilter_t *KalmanMadgwickAlloc(Coordinate_t *coor_, double timeStamp)
{
    KalmanMadgwickFilter_t *f = (KalmanMadgwickFilter_t *)malloc(sizeof(KalmanMadgwickFilter_t));
    assert(f);
    f->acc = MatrixAlloc(3, 1);
    f->gyro = MatrixAlloc(3, 1);
    f->coor = MatrixAlloc(3, 1);

    f->Rot = MatrixAlloc(3, 3); // rotation matrix
    f->E = MatrixAlloc(3, 3);   // rotation matrix
    f->kf = KalmanFilterCreate(STATE_COUNT, 3, 0);
    /*initialization*/
    f->predictTime = f->updateTime = timeStamp;

    MatrixSet(f->kf->Xk_k,
              coor_->x, coor_->y, coor_->z);

    f->kf->H->data[0][0] = 1;
    f->kf->H->data[1][1] = 1;
    f->kf->H->data[2][2] = 1;

    MatrixSetIdentity(f->kf->Pk_k);

    f->kf->Pk_k->data[0][0] = 0.5;
    f->kf->Pk_k->data[1][1] = 0.5;
    f->kf->Pk_k->data[2][2] = 0.5;

    f->kf->Pk_k->data[6][6] = 10;
    f->kf->Pk_k->data[7][7] = 10;
    f->kf->Pk_k->data[8][8] = 10;

    f->kf->Pk_k->data[9][9] = 500;
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

static void matrixAcc(KalmanMadgwickFilter_t *f, Accelerometer_t *acc_in)
{
    MatrixSet(f->acc, acc_in->x, acc_in->y, acc_in->z);
}

static void matrixGyr(KalmanMadgwickFilter_t *f, Gyroscope_t *gyr_in)
{
    MatrixSet(f->gyro, gyr_in->x, gyr_in->y, gyr_in->z);
}

static void matrixCoor(KalmanMadgwickFilter_t *f, Coordinate_t *coor_in)
{
    MatrixSet(f->coor, coor_in->x, coor_in->y, coor_in->z);
}

//////////////////////////////////////////////////////////////////////////

static void rebuildRoute(KalmanMadgwickFilter_t *f)
{
    MatrixSet(f->kf->route,
              f->kf->Xk_k->data[0][0], f->kf->Xk_k->data[1][0], f->kf->Xk_k->data[2][0]);
}

static void rebuildSpeed(KalmanMadgwickFilter_t *f)
{
    MatrixSet(f->kf->route,
              f->kf->Xk_k->data[3][0], f->kf->Xk_k->data[4][0], f->kf->Xk_k->data[5][0]);
}

static void rebuildBiasG(KalmanMadgwickFilter_t *f)
{
    MatrixSet(f->kf->bg,
              f->kf->Xk_k->data[9][0], f->kf->Xk_k->data[10][0], f->kf->Xk_k->data[11][0]);
}

static void rebuildBiasA(KalmanMadgwickFilter_t *f)
{
    MatrixSet(f->kf->ba,
              f->kf->Xk_k->data[12][0], f->kf->Xk_k->data[13][0], f->kf->Xk_k->data[14][0]);
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
    phi0 = atan2(f->Rot->data[2][1], f->Rot->data[2][2]);
    phi1 = -atan(f->Rot->data[2][0] / sqrt(1 - pow(f->Rot->data[2][0], 2)));
    phi2 = atan2(f->Rot->data[1][0], f->Rot->data[0][0]);

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

static void rebuildR(KalmanMadgwickFilter_t *f)
{
    //  MatrixSetIdentity(f->kf->R);
    //  MatrixScale(f->kf->R, posSigma*posSigma);
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

void KalmanMadgwickPredict(KalmanMadgwickFilter_t *k, double timeNow, Gyroscope_t *gyr_, Accelerometer_t *acc_)
{
    double dt = (timeNow - k->predictTime) / 1.0e+3;

    rebuildF(k, dt);
    rebuildQ(k);

    rebuildBiasA(k);
    rebuildBiasG(k);

    acc_->x -= k->kf->ba->data[0][0];
    acc_->y -= k->kf->ba->data[1][0];
    acc_->z -= k->kf->ba->data[2][0];

    gyr_->x -= k->kf->bg->data[0][0];
    gyr_->y -= k->kf->bg->data[1][0];
    gyr_->z -= k->kf->bg->data[2][0];

    matrixAcc(k, acc_);
    matrixGyr(k, gyr_);

    rebuildRoute(k);
    rebuildSpeed(k);

    MadgwickAHRSupdateIMU(gyr_->x, gyr_->y, gyr_->z, acc_->x, acc_->y, acc_->z);

    rebuildRot(k);
    rebuildE(k);

    MatrixMultiply(k->Rot, k->acc, k->kf->acc_rot);
    MatrixScale(k->kf->acc_rot, dt);
    MatrixAddIncrease(k->kf->speed, k->kf->acc_rot);

    k->kf->Xk_k->data[3][0] = k->kf->speed->data[0][0];
    k->kf->Xk_k->data[4][0] = k->kf->speed->data[1][0];
    k->kf->Xk_k->data[5][0] = k->kf->speed->data[2][0];

    MatrixScale(k->kf->speed, dt);
    MatrixMultiply(k->Rot, k->acc, k->kf->acc_rot);
    MatrixScale(k->kf->acc_rot, dt * dt * 0.5);

    MatrixAddIncrease(k->kf->route, k->kf->acc_rot);
    MatrixAddIncrease(k->kf->route, k->kf->speed);

    k->kf->Xk_k->data[0][0] = k->kf->route->data[0][0];
    k->kf->Xk_k->data[1][0] = k->kf->route->data[1][0];
    k->kf->Xk_k->data[2][0] = k->kf->route->data[2][0];

    k->predictTime = timeNow;
    KalmanFilterPredict(k->kf);
    // MatrixCopy(k->kf->Xk_km1, k->kf->Xk_k);
}
//////////////////////////////////////////////////////////////////////////

void KalmanMadgwickUpdate(KalmanMadgwickFilter_t *k, double timeNow, Coordinate_t *corr_)
{
    //  double dt = timeNow - k->updateTime;
    k->predictCount = 0;
    k->updateTime = timeNow;
    rebuildR(k);
    matrixCoor(k, corr_);
    k->kf->Zk = k->coor;

    KalmanFilterUpdate(k->kf);
}
//////////////////////////////////////////////////////////////////////////

void KalmanMadgwickGetCorr(const KalmanMadgwickFilter_t *k, Coordinate_t *corr)
{
    corr->x = k->kf->Xk_k->data[0][0];
    corr->y = k->kf->Xk_k->data[1][0];
    corr->z = k->kf->Xk_k->data[2][0];
}
//////////////////////////////////////////////////////////////////////////