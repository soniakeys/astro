// Public domain 2014, Smithsonian Astrophysical Observatory.

package astro_test

import (
	"fmt"
	"log"

	"github.com/soniakeys/astro"
)

func ExampleSolarPositionJ2000() {
	mjd := 56891.9
	e, err := astro.LoadPlanet(astro.Earth)
	if err != nil {
		log.Fatal(err)
	}
	x, y, z, _ := astro.SolarPositionJ2000(e, mjd+astro.JMod)
	fmt.Println("mjd:", mjd)
	fmt.Printf("x, y, z: %.2f %.2f %.2f\n", x, y, z)

	jd := 2448908.5
	x, y, z, r := astro.SolarPositionJ2000(e, jd)
	fmt.Printf("jd: %.1f\n", jd)
	fmt.Printf("x, y, z, r: %.6f %.6f %.6f %.6f\n", x, y, z, r)
	// Output:
	// mjd: 56891.9
	// x, y, z: -0.87 0.47 0.20
	// jd: 2448908.5
	// x, y, z, r: -0.937397 -0.313167 -0.135778 0.997609
}
