#ifndef KALMANMADGWICK_H
#define KALMANMADGWICK_H

#include <stdint.h>
#include "Kalman.h"
#include "Matrix.h"

typedef struct KalmanMadgwickFilter {
  double predictTime;
  double updateTime;
  double accDev;
  uint32_t predictCount;
  matrix_t *Rot;
  matrix_t *E;
  KalmanFilter_t *kf;
} KalmanMadgwickFilter_t;

KalmanMadgwickFilter_t* KalmanMadgwickAlloc(double x, double y,
    double xVel, double yVel,
    double accDev, double posDev,
    double timeStamp);

void KalmanMadgwickFree(KalmanMadgwickFilter_t *k);

void KalmanMadgwickPredict(KalmanMadgwickFilter_t *k, double timeNow, double xAcc, double yAcc);
void KalmanMadgwickUpdate(KalmanMadgwickFilter_t *k, double timeStamp, double x, double y, double xVel, double yVel, double posDev);

double KalmanMadgwickGetX(const KalmanMadgwickFilter_t *k);
double KalmanMadgwickGetY(const KalmanMadgwickFilter_t *k);
double KalmanMadgwickGetXVel(const KalmanMadgwickFilter_t *k);
double KalmanMadgwickGetYVel(const KalmanMadgwickFilter_t *k);

#endif // KALMANMADGWICK_H