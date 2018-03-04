// Public domain.

package astro_test

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"

	"github.com/soniakeys/astro"
	"github.com/soniakeys/meeus/v3/julian"
)

func ExampleMeeusCalendarGregorianToJD() {
	// Meeus example 7.a, p. 61.
	jd := astro.MeeusCalendarGregorianToJD(1957, 10, 4.81)
	fmt.Printf("%.2f\n", jd)
	// Meeus example 7.c, p. 64.
	jd1 := astro.MeeusCalendarGregorianToJD(1910, 4, 20)
	jd2 := astro.MeeusCalendarGregorianToJD(1986, 2, 9)
	fmt.Println(jd2 - jd1)
	// Output:
	// 2436116.31
	// 27689
}

func ExampleFFCalendarGregorianToJD() {
	// Meeus example 7.a, p. 61.
	jd := astro.FFCalendarGregorianToJD(1957, 10, 4.81)
	fmt.Printf("%.2f\n", jd)
	// Meeus example 7.c, p. 64.
	jd1 := astro.FFCalendarGregorianToJD(1910, 4, 20)
	jd2 := astro.FFCalendarGregorianToJD(1986, 2, 9)
	fmt.Println(jd2 - jd1)
	// Output:
	// 2436116.31
	// 27689
}

func ExampleTimeToJd() {
	// Meeus example 7.a, p. 61.
	t := time.Date(1957, 10, 4, 0, 0, 0, 0, time.UTC)
	ns := 0.81 * float64(24*time.Hour)
	t = t.Add(time.Duration(ns))
	fmt.Printf("%.2f\n", astro.TimeToJD(t))
	// Output:
	// 2436116.31
}

func TestFF(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	j0 := julian.CalendarGregorianToJD(1582, 11, 1)
	j1 := julian.CalendarGregorianToJD(10000, 1, 1)
	dj := j1 - j0
	for i := 0; i < 1e5; i++ {
		j := j0 + dj*rand.Float64()
		y, m, d := julian.JDToCalendar(j)
		{
			ja := astro.MeeusCalendarGregorianToJD(y, m, d)
			if dj := math.Abs(ja - j); dj > 1e-6 {
				t.Logf("%d-%02d-%f", y, m, d)
				t.Logf("ja: want %f got %f diff %f", j, ja, dj)
			}
		}
		{
			jf := astro.FFCalendarGregorianToJD(y, m, d)
			if dj := math.Abs(jf - j); dj > 1e-6 {
				t.Logf("%d-%02d-%f", y, m, d)
				t.Logf("jf: want %f got %f diff %f", j, jf, dj)
			}
		}
		{
			jm := astro.CalendarGregorianToMJD(y, m, d) + 240000.5
			if dj := math.Abs(jm - j); dj > 1e-6 {
				t.Logf("%d-%02d-%f", y, m, d)
				t.Logf("jm: want %f got %f diff %f", j, jm, dj)
			}
		}
	}
}

func BenchmarkFF(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	j0 := julian.CalendarGregorianToJD(1582, 11, 1)
	j1 := julian.CalendarGregorianToJD(2582, 11, 1)
	dj := j1 - j0
	type ymd struct {
		y, m int
		d    float64
	}
	cs := make([]ymd, 1e6)
	for i := range cs {
		y, m, d := julian.JDToCalendar(j0 + dj*rand.Float64())
		cs[i] = ymd{y, m, d}
	}
	b.Run("Meeus", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for _, c := range cs {
				astro.MeeusCalendarGregorianToJD(c.y, c.m, c.d)
			}
		}
	})
	b.Run("FF", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for _, c := range cs {
				astro.FFCalendarGregorianToJD(c.y, c.m, c.d)
			}
		}
	})
	b.Run("MJD", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for _, c := range cs {
				astro.CalendarGregorianToMJD(c.y, c.m, c.d)
			}
		}
	})
}
