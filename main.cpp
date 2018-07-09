#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <cmath>
#include "Eigen/Core"
#include "Eigen/Geometry"

using namespace Eigen;

struct Triangle
{
  Vector3d p[3];
};

struct TriangleMesh
{
  int triConut;
  std::list<Triangle> triList;
};

class Constraint
{
  public:
    int id_p1;
    int id_p2;
    double d;

    Constraint(int id1, int id2, double d)
    {
      this->id_p1 = id1;
      this->id_p2 = id2;
      this->d = d;
    }

    void calcCorrection(const std::vector<Vector3d>& pos, const std::vector<double>& invMass, Vector3d *dp1, Vector3d *dp2)
    {
      double w1 = invMass[id_p1];
      double w2 = invMass[id_p2];

      Vector3d dp12 = pos[id_p1] - pos[id_p2];
      double norm = dp12.norm();

      double delta = norm - d;
      Vector3d dir = dp12 / norm;

      (*dp1) = - w1/(w1+w2) * delta * dir;
      (*dp2) =   w2/(w1+w2) * delta * dir;
    }
};

class FloorConstraint
{
  public:
    int id;
    double y0;

    FloorConstraint(int id, double y0)
    {
      this->id = id;
      this->y0 = y0;
    }

    void calcCorrection(const std::vector<Vector3d>& pos, const std::vector<double>& invMass, Vector3d *dp)
    {
      (*dp) = Vector3d::Zero(3);

      double d = pos[id](1) - y0;
      if (d >= 0) {
        return;
      }

      (*dp)(1) = -d;
    }
};

class AreaConservationConstraint
{
  public:
    int idv[3];

    AreaConservationConstraint(int i0, int i1, int i2)
    {
      idv[0]=i0;
      idv[1]=i1;
      idv[2]=i2;
    }

    void calcCorrection(const std::vector<Vector3d>& pos, const std::vector<double>& invMass, Vector3d *dp0, Vector3d *dp1, Vector3d *dp2)
    {
      auto x10 = pos[idv[1]] - pos[idv[0]];
      auto x20 = pos[idv[2]] - pos[idv[0]];


      auto nablaX10 = 2 * x20.cross(x10.cross(x20));
      auto nablaX20 = 2 * x10.cross(x20.cross(x10));
      auto nablaX0 = - (nablaX10 + nablaX20);

      auto C = 0.25*(x10.cross(x20)).squaredNorm() - 0.25;

      auto lambda = C / (nablaX10.squaredNorm() + nablaX20.squaredNorm() + nablaX0.squaredNorm());

      (*dp1) = -lambda*nablaX10;
      (*dp2) = -lambda*nablaX20;
      (*dp0) = -lambda*nablaX0;
    }
};

class Sim {
  private:
    int N;
    std::vector<Vector3d> positions;
    std::vector<Vector3d> predictions;
    std::vector<Vector3d> velocities;
    std::vector<double> invMass;
    std::vector<Vector3d> forces;
    std::list<Constraint> constraints;
    std::list<FloorConstraint> floor_constraints;
    std::list<AreaConservationConstraint> area_constraints;

    double mass = 1.0;
    Vector3d gravity = Vector3d(0.0, -9.8, 0.0);

    double restLength = 1.0;

    double floorPosY = -7.0;

  public:
    void init(int N)
    {
      this->N = N;
      positions.resize(N);
      predictions.resize(N);
      velocities.resize(N);
      invMass.resize(N);
      forces.resize(N);
    }

    void addParticle(int id, const Vector3d &pos, const Vector3d &vel, double mass)
    {
      positions[id] = pos;
      velocities[id] = vel;
      invMass[id] = 1.0 / mass;
    }

    void addConstraint(const Constraint C)
    {
      constraints.push_back(C);
    }

    void addFloorConstraint(const FloorConstraint C)
    {
      floor_constraints.push_back(C);
    }

    void addAreaConstraint(const AreaConservationConstraint C)
    {
      area_constraints.push_back(C);
    }

    void calcForces()
    {
      // gravity
      for(int i=0; i<N; i++) {
        forces[i] = gravity / invMass[i];
      }

      /* // friction */
      /* for(int i=0; i<N; i++) { */
      /*   forces[i] += velocities[i] * (-1.0); */
      /* } */
    }

    void step(double dt)
    {
      calcForces();

      for(int i=0; i<N; i++) {
        velocities[i] += dt * invMass[i] * forces[i];
        predictions[i] = positions[i] + dt*velocities[i];
      }

      solve();

      for(int i=0; i<N; i++) {
        velocities[i] = (predictions[i] - positions[i]) / dt;
        positions[i] = predictions[i];
      }
    }

    void solve()
    {
      double eps = 0.00001 * N;
      Vector3d dp1, dp2, dp3;
      int itrMax = 100;

      double sumEps = 0;

      for(int i=0; i<itrMax; i++) {
        sumEps = 0;
        /* for(auto& c : constraints ) { */
        /*   c.calcCorrection(predictions, invMass, &dp1, &dp2); */
        /*   sumEps += dp1.squaredNorm() + dp2.squaredNorm(); */

        /*   predictions[c.id_p1] += dp1; */
        /*   predictions[c.id_p2] += dp2; */
        /* } */

        for(auto& c : area_constraints ) {
          c.calcCorrection(predictions, invMass, &dp1, &dp2, &dp3);
          sumEps += dp1.squaredNorm();

          predictions[c.idv[0]] += dp1;
          predictions[c.idv[1]] += dp2;
          predictions[c.idv[2]] += dp3;
        }

        for(auto& c : floor_constraints ) {
          c.calcCorrection(predictions, invMass, &dp1);
          sumEps += dp1.squaredNorm();

          predictions[c.id] += dp1;
        }

        if (sumEps < eps) {
          return;
        }
      }
    }

    void print()
    {
      for(int i=0; i<N; i++) {
        std::cout << positions[i](0) << " " << positions[i](1) << " " << positions[i](2) << std::endl;
      }
    }

    void write(int n)
    {
      std::ofstream ofs;
      ofs.open("dat/pos" + std::to_string(n));
      for(int i=0; i<N; i++) {
        ofs << positions[i](0) << " " << positions[i](1) << " " << positions[i](2) << std::endl;
      }
      ofs.close();
    }
};

int main(int argc, char const* argv[])
{
  Sim s;

  int N = 3;
  double mass = 1;

  s.init(N);

  double Y0 = 1.5;
  auto p0 = Vector3d(0, 0, 0);
  auto p1 = Vector3d(1, 0, 0);
  auto p2 = Vector3d(0, 1, 0);

  /* Matrix3d rot; */
  /* rot = AngleAxisd(M_PI/10, Vector3d::UnitZ()); */
  /* p0 = rot * p0; */
  /* p1 = rot * p1; */
  /* p2 = rot * p2; */
  /* p3 = rot * p3; */
  p0(1) += Y0;
  p1(1) += Y0;
  p2(1) += Y0;
  /* p3(1) += Y0; */

  s.addParticle(0, p0, Vector3d::Zero(3), mass);
  s.addParticle(1, p1, Vector3d::Zero(3), mass);
  s.addParticle(2, p2, Vector3d::Zero(3), mass);
  /* s.addParticle(3, p3, Vector3d::Zero(3), mass); */

  /* s.addConstraint(Constraint(0, 1, 1.0)); */
  /* s.addConstraint(Constraint(1, 2, sqrt(2))); */
  /* s.addConstraint(Constraint(2, 0, 1.0)); */
  /* s.addConstraint(Constraint(2, 3, 1.0)); */

  s.addFloorConstraint(FloorConstraint(0, 0));
  s.addFloorConstraint(FloorConstraint(1, 0));
  s.addFloorConstraint(FloorConstraint(2, 0));
  /* s.addFloorConstraint(FloorConstraint(3, 0)); */

  s.addAreaConstraint(AreaConservationConstraint(0, 1, 2));

  double dt = 0.01;

  s.print();

  for(int i=0; i<2000; i++) {
    s.step(dt);
    s.write(i);
  }


  return 0;
}
