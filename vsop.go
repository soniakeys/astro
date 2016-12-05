// Public domain

package astro

import (
	"errors"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/soniakeys/unit"
)

// Mercury-Neptune planet constants suitable for first argument to LoadPlanet.
const (
	Mercury = iota
	Venus
	Earth
	Mars
	Jupiter
	Saturn
	Uranus
	Neptune
	nPlanets // sad practicality
)

// parallel arrays, indexed by planet constants.
var (
	// extensions of VSOP87B files
	ext = [nPlanets]string{
		"mer", "ven", "ear", "mar", "jup", "sat", "ura", "nep"}

	// planet names as found in VSOP87B files
	b7 = [nPlanets]string{
		"MERCURY",
		"VENUS  ",
		"EARTH  ",
		"MARS   ",
		"JUPITER",
		"SATURN ",
		"URANUS ",
		"NEPTUNE",
	}
)

type abc struct {
	a, b, c float64
}

type coeff [6][]abc

// V87Planet holds VSOP87 coefficients for computing planetary
// positions in spherical coorditates.
type V87Planet struct {
	l, b, r coeff
}

// code tested with version 2.  other versions unknown.
const fileVersion = '2'

// LoadPlanet constructs a V87Planet object from a VSOP87 file.
//
// Argument ibody should be one of the planet constants.
//
// The directory containing the VSOP87 must be indicated by environment
// variable VSOP87.
func LoadPlanet(ibody int) (*V87Planet, error) {
	path := os.Getenv("VSOP87")
	if path == "" {
		return nil, errors.New("No path assigned to environment variable VSOP87")
	}
	return LoadPlanetPath(ibody, path)
}

// LoadPlanetPath constructs a V87Planet object from a VSOP87 file.
//
// Argument ibody should be one of the planet constants; path should be
// a directory containing the VSOP87 files.
func LoadPlanetPath(ibody int, path string) (*V87Planet, error) {
	if ibody < 0 || ibody >= nPlanets {
		return nil, errors.New("Invalid planet.")
	}
	data, err := ioutil.ReadFile(path + "/VSOP87B." + ext[ibody])
	if err != nil {
		return nil, err
	}
	v := &V87Planet{}
	lines := strings.Split(string(data), "\n")
	n := 0
	n, err = v.l.parse('1', ibody, lines, n, false)
	if err != nil {
		return nil, err
	}
	n, err = v.b.parse('2', ibody, lines, n, false)
	if err != nil {
		return nil, err
	}
	n, err = v.r.parse('3', ibody, lines, n, true)
	if err != nil {
		return nil, err
	}
	return v, nil
}

func (c *coeff) parse(ic byte, ibody int, lines []string, n int, au bool) (int, error) {
	var cbuf [2047]abc
	for n < len(lines) {
		line := lines[n]
		if len(line) < 132 {
			break
		}
		if line[41] != ic {
			break
		}
		if iv := line[17]; iv != fileVersion {
			return n, fmt.Errorf("Line %d: expected version %c, "+
				"found %c.", n+1, fileVersion, iv)
		}
		if bo := line[22:29]; bo != b7[ibody] {
			return n, fmt.Errorf("Line %d: expected body %s, "+
				"found %s.", n+1, b7[ibody], bo)
		}
		it := line[59] - '0'
		in, err := strconv.Atoi(strings.TrimSpace(line[60:67]))
		if err != nil {
			return n, fmt.Errorf("Line %d: %v.", n+1, err)
		}
		if in == 0 {
			continue
		}
		if in > len(lines)-n {
			return n, errors.New("Unexpected end of file.")
		}
		n++
		cx := 0
		for _, line := range lines[n : n+in] {
			a := &cbuf[cx]
			a.a, err =
				strconv.ParseFloat(strings.TrimSpace(line[79:97]), 64)
			if err != nil {
				goto parseError
			}
			a.b, err = strconv.ParseFloat(line[98:111], 64)
			if err != nil {
				goto parseError
			}
			a.c, err =
				strconv.ParseFloat(strings.TrimSpace(line[111:131]), 64)
			if err != nil {
				goto parseError
			}
			cx++
			continue
		parseError:
			return n, fmt.Errorf("Line %d: %v.", n+cx+1, err)
		}
		c[it] = append([]abc{}, cbuf[:cx]...)
		n += in
	}
	return n, nil
}

// Position2000 returns ecliptic position of planets by full VSOP87 theory.
//
// Argument jde is the date for which positions are desired.
//
// Results are for the dynamical equinox and ecliptic J2000.
//
//	L is heliocentric longitude in radians.
//	B is heliocentric latitude in radians.
//	R is heliocentric range in AU.
func (vt *V87Planet) Position2000(jde float64) (L, B unit.Angle, R float64) {
	T := J2000Century(jde)
	τ := T * .1
	cf := make([]float64, 6)
	sum := func(series coeff) float64 {
		for x, terms := range series {
			cf[x] = 0
			// sum terms in reverse order to preserve accuracy
			for y := len(terms) - 1; y >= 0; y-- {
				term := &terms[y]
				cf[x] += term.a * math.Cos(term.b+term.c*τ)
			}
		}
		return Horner(τ, cf[:len(series)]...)
	}
	L = unit.Angle(sum(vt.l)).Mod1()
	B = unit.Angle(sum(vt.b))
	R = sum(vt.r)
	return
}
