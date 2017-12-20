// Public domain

package astro

import "time"

// MeeusCalendarGregorianToJD converts a Gregorian year, month, and day of
// month to Julian day.
//
// Negative years are valid, back to JD 0.  The result is not valid for
// dates before JD 0.
//
// The algorithm by Meeus relies on floor truncation like the Fortran INT
// function.  The implemenation here retains the floor truncation, which
// makes this algorithm relatively slow.
func MeeusCalendarGregorianToJD(y, m int, d float64) float64 {
	switch m {
	case 1, 2:
		y--
		m += 12
	}
	a := floorDiv(y, 100)
	b := 2 - a + floorDiv(a, 4)
	// (7.1) p. 61
	return float64(floorDiv64(36525*(int64(y+4716)), 100)) +
		float64(floorDiv(306*(m+1), 10)+b) + d - 1524.5
}

// FFCalendarGregorianToJD converts a Gregorian year, month, and day of
// month to Julian day.
//
// Implementation is by the widely used formula of Fliegel and van Flandern.
// It is fast, although the formula appears rather cryptic.
func FFCalendarGregorianToJD(y, m int, d float64) float64 {
	return d - .5 + float64(-32075+1461*(y+4800+(m-14)/12)/4+
		367*(m-2-(m-14)/12*12)/12-3*((y+4900+(m-14)/12)/100)/4)
}

// CalendarGregorianToMJD converts a Gregorian year, month, and day of
// month to Modified Julian Day (MJD).
//
// The algorithm here is different than either the Meeus or Fliegel -
// van Flandern Julian date algorithms.  It is more intuitive than either
// yet runs about as fast as the Fliegel - van Flandern algorithm.
//
// It is not valid for proleptic dates with negative years.
func CalendarGregorianToMJD(y, m int, d float64) float64 {
	if m <= 2 {
		y--
	}
	return float64(int(daysSinceFeb[m])+365*y+y/4-y/100+y/400-678882) + d
}

type dsf int

var daysSinceFeb = [...]dsf{
	0, jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec}

const (
	// at the start of mar 1, there have been 0 days
	mar dsf = 0
	// for apr, add the 31 days of mar
	apr = mar + 31
	may = apr + 30
	jun = may + 31
	jul = jun + 30
	aug = jul + 31
	sep = aug + 31
	oct = sep + 30
	nov = oct + 31
	dec = nov + 30
	jan = dec + 31
	feb = jan + 31
)

// TimeToJD takes a Go time.Time and returns a JD as float64.
//
// Any time zone offset in the time.Time is ignored and the time is
// treated as UTC.
func TimeToJD(t time.Time) float64 {
	ut := t.UTC()
	y, m, _ := ut.Date()
	d := ut.Sub(time.Date(y, m, 0, 0, 0, 0, 0, time.UTC))
	// time.Time is always Gregorian
	return FFCalendarGregorianToJD(y, int(m), float64(d)/float64(24*time.Hour))
}
