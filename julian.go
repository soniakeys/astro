// Public domain

package astro

import "time"

// CalendarGregorianToJD converts a Gregorian year, month, and day of month
// to Julian day.
//
// Negative years are valid, back to JD 0.  The result is not valid for
// dates before JD 0.
func CalendarGregorianToJD(y, m int, d float64) float64 {
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

// TimeToJD takes a Go time.Time and returns a JD as float64.
//
// Any time zone offset in the time.Time is ignored and the time is
// treated as UTC.
func TimeToJD(t time.Time) float64 {
	ut := t.UTC()
	y, m, _ := ut.Date()
	d := ut.Sub(time.Date(y, m, 0, 0, 0, 0, 0, time.UTC))
	// time.Time is always Gregorian
	return CalendarGregorianToJD(y, int(m), float64(d)/float64(24*time.Hour))
}
