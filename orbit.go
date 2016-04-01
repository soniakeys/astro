// Public domain

package astro

import (
	"errors"
	"math"
)

type Elements struct {
	Axis  float64 // Semimajor axis, a, in AU
	Ecc   float64 // Eccentricity, e
	Inc   float64 // Inclination, i, in radians
	ArgP  float64 // Argument of perihelion, ω, in radians
	Node  float64 // Longitude of ascending node, Ω, in radians
	TimeP float64 // Time of perihelion, T, as jde
}
type Orbit struct {
	k          *Elements
	n          float64
	_A, _B, _C float64
	a, b, c    float64
}

func NewOrbit(k *Elements) *Orbit {
	o := &Orbit{
		k: k,
		n: K / k.Axis / math.Sqrt(k.Axis),
	}
	const sε = SOblJ2000
	const cε = COblJ2000
	sΩ, cΩ := math.Sincos(k.Node)
	si, ci := math.Sincos(k.Inc)
	// (33.7) p. 228
	F := cΩ
	G := sΩ * cε
	H := sΩ * sε
	P := -sΩ * ci
	Q := cΩ*ci*cε - si*sε
	R := cΩ*ci*sε + si*cε
	// (33.8) p. 229
	o._A = math.Atan2(F, P)
	o._B = math.Atan2(G, Q)
	o._C = math.Atan2(H, R)
	o.a = math.Hypot(F, P)
	o.b = math.Hypot(G, Q)
	o.c = math.Hypot(H, R)
	return o
}

func (o *Orbit) Position(jde float64) (x, y, z, r float64) {
	M := o.n * (jde - o.k.TimeP)
	E := kepler(o.k.Ecc, M)
	r = radius(E, o.k.Ecc, o.k.Axis)
	ν := trueAnomaly(E, o.k.Ecc)
	// (33.9) p. 229
	x = r * o.a * math.Sin(o._A+o.k.ArgP+ν)
	y = r * o.b * math.Sin(o._B+o.k.ArgP+ν)
	z = r * o.c * math.Sin(o._C+o.k.ArgP+ν)
	return
}

func kepler(e, M float64) float64 {
	if E, err := kepler2b(e, M, 15); err != nil {
		return E
	}
	return kepler3(e, M)
}

// Kepler2b solves Kepler's equation by iteration.
//
// The iterated formula is the same as in Kepler2 but a (different) limiting
// function avoids divergence.
//
// Argument e is eccentricity, M is mean anomaly in radians,
// places is the desired number of decimal places in the result.
//
// Result E is eccentric anomaly in radians.
func kepler2b(e, M float64, places int) (E float64, err error) {
	f := func(E0 float64) float64 {
		se, ce := math.Sincos(E0)
		d := (M + e*se - E0) / (1 - e*ce)
		// method of Steele, p. 205
		if d > .5 {
			d = .5
		} else if d < -.5 {
			d = -.5
		}
		return E0 + d
	}
	return iterateDecimalPlaces(f, M, places, places)
}

// Kepler3 solves Kepler's equation by binary search.
//
// Argument e is eccentricity, M is mean anomaly in radians.
//
// Result E is eccentric anomaly in radians.
func kepler3(e, M float64) (E float64) {
	// adapted from BASIC, p. 206
	M = PMod(M, 2*math.Pi)
	f := 1
	if M > math.Pi {
		f = -1
		M = 2*math.Pi - M
	}
	E0 := math.Pi * .5
	d := math.Pi * .25
	for i := 0; i < 53; i++ {
		M1 := E0 - e*math.Sin(E0)
		if M-M1 < 0 {
			E0 -= d
		} else {
			E0 += d
		}
		d *= .5
	}
	if f < 0 {
		return -E0
	}
	return E0
}

// DecimalPlaces iterates to a fixed number of decimal places.
//
// Inputs are an improvement function, a starting value, the number of
// decimal places desired in the result, and an iteration limit.
func iterateDecimalPlaces(better func(float64) float64, start float64, places, maxIterations int) (float64, error) {
	d := math.Pow(10, float64(-places))
	for i := 0; i < maxIterations; i++ {
		n := better(start)
		if math.Abs(n-start) < d {
			return n, nil
		}
		start = n
	}
	return 0, errors.New("Maximum iterations reached")
}

// True returns true anomaly ν for given eccentric anomaly E.
//
// Argument e is eccentricity.  E must be in radians.
//
// Result is in radians.
func trueAnomaly(E, e float64) float64 {
	// (30.1) p. 195
	return 2 * math.Atan(math.Sqrt((1+e)/(1-e))*math.Tan(E*.5))
}

// Radius returns radius distance r for given eccentric anomaly E.
//
// Argument e is eccentricity, a is semimajor axis.
//
// Result unit is the unit of semimajor axis a (typically AU.)
func radius(E, e, a float64) float64 {
	// (30.2) p. 195
	return a * (1 - e*math.Cos(E))
}
