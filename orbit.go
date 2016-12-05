// Public domain

package astro

import (
	"errors"
	"math"

	"github.com/soniakeys/unit"
)

type Elements struct {
	Axis  float64    // Semimajor axis, a, in AU
	Ecc   float64    // Eccentricity, e
	Inc   unit.Angle // Inclination, i
	ArgP  unit.Angle // Argument of perihelion, ω
	Node  unit.Angle // Longitude of ascending node, Ω
	TimeP float64    // Time of perihelion, T, as jde
}
type Orbit struct {
	k          *Elements
	n          unit.Angle // Angle/day
	_A, _B, _C unit.Angle
	a, b, c    float64
}

func NewOrbit(k *Elements) *Orbit {
	o := &Orbit{
		k: k,
		n: unit.Angle(K / k.Axis / math.Sqrt(k.Axis)),
	}
	const sε = SOblJ2000
	const cε = COblJ2000
	sΩ, cΩ := k.Node.Sincos()
	si, ci := k.Inc.Sincos()
	// (33.7) p. 228
	F := cΩ
	G := sΩ * cε
	H := sΩ * sε
	P := -sΩ * ci
	Q := cΩ*ci*cε - si*sε
	R := cΩ*ci*sε + si*cε
	// (33.8) p. 229
	o._A = unit.Angle(math.Atan2(F, P))
	o._B = unit.Angle(math.Atan2(G, Q))
	o._C = unit.Angle(math.Atan2(H, R))
	o.a = math.Hypot(F, P)
	o.b = math.Hypot(G, Q)
	o.c = math.Hypot(H, R)
	return o
}

func (o *Orbit) Position(jde float64) (x, y, z, r float64) {
	M := o.n.Mul(jde - o.k.TimeP)
	E := kepler(o.k.Ecc, M)
	r = radius(E, o.k.Ecc, o.k.Axis)
	ν := trueAnomaly(E, o.k.Ecc)
	// (33.9) p. 229
	x = r * o.a * (o._A + o.k.ArgP + ν).Sin()
	y = r * o.b * (o._B + o.k.ArgP + ν).Sin()
	z = r * o.c * (o._C + o.k.ArgP + ν).Sin()
	return
}

func kepler(e float64, M unit.Angle) unit.Angle {
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
func kepler2b(e float64, M unit.Angle, places int) (E unit.Angle, err error) {
	f := func(E0 float64) float64 {
		se, ce := math.Sincos(E0)
		d := (M.Rad() + e*se - E0) / (1 - e*ce)
		// method of Steele, Meeus p. 205
		if d > .5 {
			d = .5
		} else if d < -.5 {
			d = -.5
		}
		return E0 + d
	}
	ea, err := iterateDecimalPlaces(f, M.Rad(), places, places)
	return unit.Angle(ea), err
}

// Kepler3 solves Kepler's equation by binary search.
//
// Argument e is eccentricity, M is mean anomaly in radians.
//
// Result E is eccentric anomaly in radians.
func kepler3(e float64, M unit.Angle) (E unit.Angle) {
	// adapted from BASIC, Meeus p. 206
	MR := M.Mod1().Rad()
	f := 1
	if MR > math.Pi {
		f = -1
		MR = 2*math.Pi - MR
	}
	E0 := math.Pi * .5
	d := math.Pi * .25
	for i := 0; i < 53; i++ {
		M1 := E0 - e*math.Sin(E0)
		if MR-M1 < 0 {
			E0 -= d
		} else {
			E0 += d
		}
		d *= .5
	}
	if f < 0 {
		E0 = -E0
	}
	return unit.Angle(E0)
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
// Argument e is eccentricity.
func trueAnomaly(E unit.Angle, e float64) unit.Angle {
	// (30.1) p. 195
	return unit.Angle(2 * math.Atan(math.Sqrt((1+e)/(1-e))*E.Mul(.5).Tan()))
}

// Radius returns radius distance r for given eccentric anomaly E.
//
// Argument e is eccentricity, a is semimajor axis.
//
// Result unit is the unit of semimajor axis a (typically AU.)
func radius(E unit.Angle, e, a float64) float64 {
	// (30.2) p. 195
	return a * (1 - e*E.Cos())
}
