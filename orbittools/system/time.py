import numpy as np


def julian_date(year, month, day, hour, minute, second):
    return (367 * year) - int((year + int((month+9)/12))*7/4) + int(275*month/9) + day + 1721013.5 + (((((second/60)+minute)/60)+hour)/24)


def gregorian_date(julian_date):
    T1900 = (julian_date - 2415019.5) / 365.25
    year = 1900 + np.floor(T1900)
    leap_years = np.floor((year - 1900 - 1) / 4)
    days = (julian_date - 2415019.5) - (((year - 1900) * 365) + leap_years)
    if days < 1:
        year = year - 1
        leap_years = np.floor((year - 1900 - 1) / 4)
        days = (julian_date - 2415019.5) - (((year - 1900) * 365) + leap_years)

    month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if year % 4 == 0:
        febdays = 29
    day_of_year = np.floor(days)
    sum_of_days = 0
    month = 1
    for i, month_days in enumerate(month_lengths):
        sum_of_days += month_days
        month += 1
        if sum_of_days + 1 > day_of_year:
            break
    day = day_of_year - sum_of_days
    time = 24 * (days - day_of_year)
    hour = np.floor(time)
    minute = np.floor((time - hour) * 60)
    second = (time - hour - (minute/60)) * 3600
    return year, month, day, hour, minute, second
