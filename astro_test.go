// Public domain 2014, Smithsonian Astrophysical Observatory.

package astro_test

import (
	"fmt"
	"math"

	"github.com/soniakeys/astro"
	"github.com/soniakeys/coord"
	"github.com/soniakeys/sexagesimal"
)

func ExampleAeiHv() {
	fmt.Println("Example 1")
	p := &coord.Cart{1, 0, 0}
	v := &coord.Cart{0, 1, 0}
	f := func() {
		a, e, i, hpv := astro.AeiHv(p, v, math.Sqrt(p.Square()))
		fmt.Printf("a: %.5f AU\n", a)
		fmt.Printf("e: %.5f\n", e)
		fmt.Printf("i: %.5f deg\n", i)
		fmt.Printf("hpv: {X:%.2f Y:%.2f Z:%.2f}\n", hpv.X, hpv.Y, hpv.Z)
	}
	f()

	fmt.Print("\nExample 2\n")
	p = &coord.Cart{1.5, 1.5, .2}
	v = &coord.Cart{-.5, .5, 0}
	f()
	// Output:
	// Example 1
	// a: 1.00000 AU
	// e: 0.00000
	// i: 0.00000 deg
	// hpv: {X:0.00 Y:0.00 Z:1.00}
	//
	// Example 2
	// a: 2.27974 AU
	// e: 0.06536
	// i: 5.38598 deg
	// hpv: {X:-0.10 Y:-0.10 Z:1.50}
}

func ExampleHMag() {
	oov := &coord.Cart{.5, 0, 0}
	sov := &coord.Cart{1, 0, 0}
	vmag := 20.
	ood := .5
	sod := 1.
	fmt.Printf("H = %.1f\n", astro.HMag(oov, sov, vmag, ood, sod))
	// Output:
	// H = 21.5
}

func ExampleLst() {
	mjd := 46895.
	fmt.Println(sexa.NewFmtHourAngle(astro.Lst(mjd, 0)))
	fmt.Println(sexa.NewFmtHourAngle(astro.Lst(mjd, .01)))
	// Output:
	// 13ʰ10ᵐ46ˢ
	// 13ʰ25ᵐ10ˢ
}

// ExampleLst output in rough agreement with Meeus, 2nd ed p 88,
// or implementation github.com/soniakeys/meeus/sidereal.Mean().

func ExampleSe2000() {
	mjd := 56891.9
	sunEarth, s, c := astro.Se2000(mjd)
	fmt.Printf("{X:%.3f Y:%.3f Z:%.3f}\n",
		sunEarth.X, sunEarth.Y, sunEarth.Z)
	fmt.Printf("%.1f\n", math.Atan2(s, c)*180/math.Pi)
	// Output:
	// {X:-0.873 Y:0.468 Z:0.203}
	// 23.4
}

// ExampleSe2000 output in rough agreement with aprx.Equ(EclPos)

func ExamplePMod() {
	fmt.Println(math.Mod(-1, 3))
	fmt.Println(astro.PMod(-1, 3))
	fmt.Println(astro.PMod(-1, -3))
	// Output:
	// -1
	// 2
	// -4
}
