// Public domain.

// Package astro, stuff generally useful in astronomy.
package astro

import (
	"math"

	"github.com/soniakeys/coord"
)

const (
	C    = 299792458    // Speed of light, m/s
	AU   = 149597870691 // 1 AU in meters
	K    = .01720209895 // Gaussian gravitational constant
	InvK = 1 / K
	U    = K * K
)

// SOblJ2000, COblJ2000 are sine and cosine of obliquity at J2000.
const (
	SOblJ2000 = .397777156
	COblJ2000 = .917482062
)

const (
	J2000         = 2451545.0 // Julian date corresponding to January 1.5, year 2000.
	JulianCentury = 36525     // days
)

// J2000Century returns the number of Julian centuries since J2000.
//
// The quantity appears as T in a number of time series.
func J2000Century(jde float64) float64 {
	return (jde - J2000) / JulianCentury
}

// AeiHv, solves Keplerian elements from state vectors.
//
// Actually stretching the package claim of "generally useful," this
// function has parameters and return values most efficient for Digest2.
//
// The algorithm becomes unstable for near-parabolic orbits or orbits with
// large semimajor axes.  The function returns ok=false if a would be > 100 AU
// or if e would be > .99.
//
// Args:
//   p = position: sun object vector, in AU
//   v = object velocity vector, scaled by gravitational constant
//   d = sun-object distance pre computed from p
//   hv: a result argument, must be a non-nil pointer, used to return
//       the momentum vector.
//
// Returns keplerian a, e, and i, the momentum vector (through an argument),
// and ok=true for valid results.
//
// Keplerian a returned in AU, i returned in degrees.
//
// ok=false means other return values are invalid
func AeiHv(p, v *coord.Cart, d float64, hv *coord.Cart) (a, e, i float64, ok bool) {

	// momentum vector
	hv.Cross(p, v)
	hsq := hv.Square()
	hm := math.Sqrt(hsq)

	// solve for semi-major axis
	// (and the inverse--it comes in handy)
	vsq := v.Square()
	temp := 2 - d*vsq

	// for stability, require a < 100
	if d > temp*100 {
		return // with ok=false
	}
	a = d / temp
	inva := temp / d

	// solve for eccentricity
	// (stability test on a (above) should keep result real)
	e = math.Sqrt(1 - hsq*inva)

	// stability test:  require e < .99
	if e > .99 {
		return // with ok=false
	}

	// solve for inclination.

	// reliable check for i=0.  handles loss of precision in h computation.
	iZero := hv.Z >= hm
	// combination of stability tests on a and e (above) should
	// ensure that hm is well above zero.
	if !iZero {
		i = math.Acos(hv.Z/hm) * 180 / math.Pi
	}
	return a, e, i, true
}

// HMag computes H from V magnitude.
func HMag(oov, sov *coord.Cart, vmag, ood, sod float64) float64 {
	rdelta := ood * sod
	cospsi := oov.Dot(sov) / rdelta

	if cospsi < -.9999 {
		// object is straight into the sun.  doesn't seem too likely,
		// but anyway, this returns a valid value.
		return 30
	}

	tanhalf := math.Sqrt(1-cospsi*cospsi) / (1 + cospsi)
	phi1 := math.Exp(-3.33 * math.Pow(tanhalf, 0.63))
	phi2 := math.Exp(-1.87 * math.Pow(tanhalf, 1.22))
	return vmag -
		5.*math.Log10(rdelta) +
		2.5*math.Log10(.85*phi1+.15*phi2)
}

// Lst computes (approximate) local sidereal time.
//
// Argument mjd is modified Julian day, long is longitude in circles.
//
//  Returns local sidereal time in radians where 1 day = 2π radians.
func Lst(mjd, long float64) float64 {
	t := (mjd - 15019.5) / 36525
	th := (6.6460656 + (2400.051262+0.00002581*t)*t) / 24
	_, ut := math.Modf(mjd)
	if ut < 0 {
		ut++
	}
	_, s := math.Modf(th + ut + long)
	if s < 0 {
		s++
	}
	return s * 2 * math.Pi
}

// Se2000 computes solar ephemeris, J2000.
//
// Returns:
//   sunEarth:  sun-earth vector in equatorial coordinates.
//   soe, coe:  sine and cosine of ecciptic.
//
// Notes:
//   Approximate solar coordinates, per USNO.  Originally from
//   http://aa.usno.navy.mil/faq/docs/SunApprox.html, page now at
//   http://www.usno.navy.mil/USNO/astronomical-applications/
//   astronomical-information-center/approx-solar.
func Se2000(mjd float64) (sunEarth coord.Cart, soe, coe float64) {
	// USNO algorithm is in degrees.  To mimimize confusion, work in
	// degrees here too, only converting to radians as needed for trig
	// functions.
	d := mjd - 51544.5
	g := 357.529 + .98560028*d // mean anomaly of sun, in degrees
	q := 280.459 + .98564736*d // mean longitude of sun, in degrees
	g2 := g + g
	sg, cg := math.Sincos(g * math.Pi / 180) // send radians to trig function
	sg2, cg2 := math.Sincos(g2 * math.Pi / 180)

	// ecliptic longitude, in degrees still
	l := q + 1.915*sg + .020*sg2

	// distance in AU
	r := 1.00014 - .01671*cg - .00014*cg2

	// obliquity of ecliptic in degrees
	e := 23.439 - .00000036*d
	soe, coe = math.Sincos(e * math.Pi / 180)

	// equatorial coordinates
	sl, cl := math.Sincos(l * math.Pi / 180)
	sunEarth.X = r * cl
	rsl := r * sl
	sunEarth.Y = rsl * coe
	sunEarth.Z = rsl * soe
	return
}

// Horner evaluates a polynomal with coefficients c at x.  The constant
// term is c[0].  The function panics with an empty coefficient list.
func Horner(x float64, c ...float64) float64 {
	i := len(c) - 1
	y := c[i]
	for i > 0 {
		i--
		y = y*x + c[i] // sorry, no fused multiply-add in Go
	}
	return y
}

// PMod returns a positive floating-point x mod y.
//
// For a positive argument y, it returns a value in the range [0,y).
//
// The function is not useful for negative values of y.
func PMod(x, y float64) float64 {
	r := math.Mod(x, y)
	if r < 0 {
		r += y
	}
	return r
}

// FloorDiv returns the integer floor of the fractional value (x / y).
//
// It uses integer math only, so is more efficient than using floating point
// intermediate values.  This function can be used in many places where INT()
// appears in AA.  As with built in integer division, it panics with y == 0.
func floorDiv(x, y int) int {
	if (x < 0) == (y < 0) {
		return x / y
	}
	return x/y - 1
}

// FloorDiv64 returns the integer floor of the fractional value (x / y).
//
// It uses integer math only, so is more efficient than using floating point
// intermediate values.  This function can be used in many places where INT()
// appears in AA.  As with built in integer division, it panics with y == 0.
func floorDiv64(x, y int64) int64 {
	if (x < 0) == (y < 0) {
		return x / y
	}
	return x/y - 1
}
