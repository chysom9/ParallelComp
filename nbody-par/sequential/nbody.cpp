#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
#include <string>
#include <omp.h>
// #include <omp.h>

double G = 6.674 * std::pow(10, -11);
double SOFTENING = 1e-1;

// double G = 1;

struct simulation
{
  size_t nbpart;

  std::vector<double> mass;

  // position
  std::vector<double> x, y, z;

  // velocity
  std::vector<double> vx, vy, vz;

  // force
  std::vector<double> fx, fy, fz;

  simulation(size_t nb)
      : nbpart(nb), mass(nb),
        x(nb), y(nb), z(nb),
        vx(nb), vy(nb), vz(nb),
        fx(nb), fy(nb), fz(nb)
  {
  }
};

void random_init(simulation &s)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution dismass(0.9, 1.);
  std::normal_distribution dispos(0., 1.);
  std::normal_distribution disvel(0., 1.);

  for (size_t i = 0; i < s.nbpart; ++i)
  {
    s.mass[i] = dismass(gen);

    s.x[i] = dispos(gen);
    s.y[i] = dispos(gen);
    s.z[i] = dispos(gen); // confused why this is here but ill leave it in
    s.z[i] = 0.;

    s.vx[i] = disvel(gen); // confused why this is here but ill leave it in
    s.vy[i] = disvel(gen); // confused why this is here but ill leave it in
    s.vz[i] = disvel(gen); // confused why this is here but ill leave it in
    s.vz[i] = 0.;
    s.vx[i] = s.y[i] * 1.5;
    s.vy[i] = -s.x[i] * 1.5;
  }

  return; //???

  // normalize velocity (using normalization found on some physicis blog)
  double meanmass = 0;
  double meanmassvx = 0;
  double meanmassvy = 0;
  double meanmassvz = 0;
  for (size_t i = 0; i < s.nbpart; ++i)
  {
    meanmass += s.mass[i];
    meanmassvx += s.mass[i] * s.vx[i];
    meanmassvy += s.mass[i] * s.vy[i];
    meanmassvz += s.mass[i] * s.vz[i];
  }
  for (size_t i = 0; i < s.nbpart; ++i)
  {
    s.vx[i] -= meanmassvx / meanmass;
    s.vy[i] -= meanmassvy / meanmass;
    s.vz[i] -= meanmassvz / meanmass;
  }
}

void init_solar(simulation &s)
{
  enum Planets
  {
    SUN,
    MERCURY,
    VENUS,
    EARTH,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    MOON
  };
  s = simulation(10);

  // Masses in kg
  s.mass[SUN] = 1.9891 * std::pow(10, 30);
  s.mass[MERCURY] = 3.285 * std::pow(10, 23);
  s.mass[VENUS] = 4.867 * std::pow(10, 24);
  s.mass[EARTH] = 5.972 * std::pow(10, 24);
  s.mass[MARS] = 6.39 * std::pow(10, 23);
  s.mass[JUPITER] = 1.898 * std::pow(10, 27);
  s.mass[SATURN] = 5.683 * std::pow(10, 26);
  s.mass[URANUS] = 8.681 * std::pow(10, 25);
  s.mass[NEPTUNE] = 1.024 * std::pow(10, 26);
  s.mass[MOON] = 7.342 * std::pow(10, 22);

  // Positions (in meters) and velocities (in m/s)
  double AU = 1.496 * std::pow(10, 11); // Astronomical Unit

  s.x = {0, 0.39 * AU, 0.72 * AU, 1.0 * AU, 1.52 * AU, 5.20 * AU, 9.58 * AU, 19.22 * AU, 30.05 * AU, 1.0 * AU + 3.844 * std::pow(10, 8)};
  s.y = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  s.z = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  s.vx = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  s.vy = {0, 47870, 35020, 29780, 24130, 13070, 9680, 6800, 5430, 29780 + 1022};
  s.vz = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
}

// meant to update the force that from applies on to
// deleted update force because that is actually confusing to read lol
//  could probably be written more concisely.

void reset_force(simulation &s)
{
#pragma omp parallel for
  for (size_t i = 0; i < s.nbpart; ++i)
  {
    s.fx[i] = 0.;
    s.fy[i] = 0.;
    s.fz[i] = 0.;
  }
}

void compute_forces(simulation &s)
{
// run this whole block in parallel threads
#pragma omp parallel
  {
    // Each thread gets its own temporary force arrays
    // so they don’t step on each other
    std::vector<double> fx_local(s.nbpart, 0.0);
    std::vector<double> fy_local(s.nbpart, 0.0);
    std::vector<double> fz_local(s.nbpart, 0.0);

// Split the outer loop (over each particle i) across threads
#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < s.nbpart; ++i)
    {
      // For particle i, look at every other particle j
      for (size_t j = 0; j < s.nbpart; ++j)
      {
        if (i == j)
          continue; // skip self

        // 1) Compute the vector from i to j
        double dx = s.x[j] - s.x[i];
        double dy = s.y[j] - s.y[i];
        double dz = s.z[j] - s.z[i];

        // 2) Compute squared distance + a little softening
        double dist2 = dx * dx + dy * dy + dz * dz + SOFTENING;

        // 3) Inverse distance and inverse-distance-squared
        double invDist = 1.0 / std::sqrt(dist2);
        double invDist3 = invDist * invDist * invDist;

        // 4) Gravitational force magnitude
        //    (G * m_i * m_j) / r^2  becomes  G * m_i * m_j * invDist^2
        double F = G * s.mass[i] * s.mass[j] * (invDist * invDist);

        // 5) Accumulate force components onto our local buffer
        fx_local[i] += F * dx * invDist;
        fy_local[i] += F * dy * invDist;
        fz_local[i] += F * dz * invDist;
      }
    }

// Now merge each thread’s local forces into the global arrays
// We do this one particle at a time under a lock
#pragma omp critical
    {
      for (size_t i = 0; i < s.nbpart; ++i)
      {
        s.fx[i] += fx_local[i];
        s.fy[i] += fy_local[i];
        s.fz[i] += fz_local[i];
      }
    }
  }
}

void integrate(simulation &s, double dt)
{
#pragma omp parallel for
  for (size_t i = 0; i < s.nbpart; ++i)
  {
    s.vx[i] += (s.fx[i] / s.mass[i]) * dt;
    s.vy[i] += (s.fy[i] / s.mass[i]) * dt;
    s.vz[i] += (s.fz[i] / s.mass[i]) * dt;
    s.x[i] += s.vx[i] * dt;
    s.y[i] += s.vy[i] * dt;
    s.z[i] += s.vz[i] * dt;
  }
}

void dump_state(simulation &s)
{
  std::cout << s.nbpart << '\t';
  for (size_t i = 0; i < s.nbpart; ++i)
  {
    std::cout << s.mass[i] << '\t';
    std::cout << s.x[i] << '\t' << s.y[i] << '\t' << s.z[i] << '\t';
    std::cout << s.vx[i] << '\t' << s.vy[i] << '\t' << s.vz[i] << '\t';
    std::cout << s.fx[i] << '\t' << s.fy[i] << '\t' << s.fz[i] << '\t';
  }
  std::cout << '\n';
}

void load_from_file(simulation &s, std::string filename)
{
  std::ifstream in(filename);
  size_t nbpart;
  in >> nbpart;
  s = simulation(nbpart);
  for (size_t i = 0; i < s.nbpart; ++i)
  {
    in >> s.mass[i];
    in >> s.x[i] >> s.y[i] >> s.z[i];
    in >> s.vx[i] >> s.vy[i] >> s.vz[i];
    in >> s.fx[i] >> s.fy[i] >> s.fz[i];
  }
  if (!in.good())
    throw "NOPEEEE"; // lol
}

int main(int argc, char *argv[])
{
  if (argc != 5)
  {
    std::cerr
        << "usage: " << argv[0] << " <input> <dt> <nbstep> <printevery>" << "\n"
        << "input can be:" << "\n"
        << "a number (random initialization)" << "\n"
        << "planet (initialize with solar system)" << "\n"
        << "a filename (load from file in singleline tsv)" << "\n";
    return -1;
  }
  std::string input(argv[1]);
  double dt = std::atof(argv[2]); // in seconds
  size_t nbstep = std::atol(argv[3]);
  size_t printevery = std::atol(argv[4]);

  simulation s(1);

  // parse command line
  // {
  //   size_t nbpart = std::atol(argv[1]); // return 0 if not a number
  //   if (nbpart > 0)
  //   {
  //     s = simulation(nbpart);
  //     random_init(s);
  //   }
  //   else
  //   {
  //     std::string inputparam = argv[1];
  //     if (inputparam == "planet")
  //     {
  //       init_solar(s);
  //     }
  //     else
  //     {
  //       load_from_file(s, inputparam);
  //     }
  //   }
  // }

  try {
    size_t n = std::stoul(input);
    s = simulation(n);
    random_init(s);
  } catch (...) {
    if (input == "planet")
      init_solar(s);
    else
      load_from_file(s, input);
  }

  for (size_t step = 0; step < nbstep; ++step) {
    if (step % printevery == 0)
      dump_state(s);

    reset_force(s);
    compute_forces(s);
    integrate(s, dt);
  }
}

