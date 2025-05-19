#ifndef KALMANMADGWICK_H
#define KALMANMADGWICK_H

#include <stdint.h>
#include "Kalman.h"
#include "Matrix.h"

typedef struct KalmanMadgwickFilter
{
  double predictTime;
  double updateTime;
  uint32_t predictCount;
  matrix_t *coor;
  matrix_t *acc;
  matrix_t *gyro;
  matrix_t *Rot;
  matrix_t *E;
  KalmanFilter_t *kf;
} KalmanMadgwickFilter_t;

KalmanMadgwickFilter_t *KalmanMadgwickAlloc(Coordinate_t *,double timeStamp);

void KalmanMadgwickFree(KalmanMadgwickFilter_t *k);

void KalmanMadgwickPredict(KalmanMadgwickFilter_t *k, double timeNow, Gyroscope_t *, Accelerometer_t *);
void KalmanMadgwickUpdate(KalmanMadgwickFilter_t *k, double timeNow, Coordinate_t *);

void KalmanMadgwickGetCorr(const KalmanMadgwickFilter_t *k, Coordinate_t *corr);

#endif // KALMANMADGWICK_H