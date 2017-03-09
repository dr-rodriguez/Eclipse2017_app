from flask import Flask, render_template, request, redirect

from astropy.time import Time
from astropy.coordinates import EarthLocation, get_sun, AltAz
from moon import moon_distance, get_umbra, estimate_umbra

app = Flask(__name__)

app.vars = dict()
app.vars['location'] = 'Nashville, TN'
app.vars['c2time'] = '18:27:03.9'
app.vars['c3time'] = '18:29:43.0'
app.vars['d_est'] = ''
app.vars['d_real'] = ''
app.vars['su_est'] = ''
app.vars['su_real'] = ''
app.vars['lat'] = ''
app.vars['long'] = ''

# Redirect to the main page
@app.route('/')
@app.route('/index.html')
@app.route('/query.html')
@app.route('/query')
def bdnyc_home():
    return redirect('/index')


# Page with a text box to take the user parameters
@app.route('/index', methods=['GET', 'POST'])
def bdnyc_query():
    return render_template('index.html', location=app.vars['location'],
                           c2time=app.vars['c2time'], c3time=app.vars['c3time'])


# Grab results of query and calculate the distance
@app.route('/runcalc', methods=['POST'])
def bdnyc_runquery():
    app.vars['location'] = request.form['location']
    app.vars['c2time'] = request.form['c2time']
    app.vars['c3time'] = request.form['c3time']

    loc = EarthLocation.of_address(app.vars['location'])
    t2 = Time('2017-08-21 ' + app.vars['c2time'])
    t3 = Time('2017-08-21 ' + app.vars['c3time'])
    dm, reald = moon_distance(loc, t2, t3)

    sun = get_sun(t2)
    sun_altaz = sun.transform_to(AltAz(obstime=t2, location=loc))

    su_est = estimate_umbra(loc, t2, t3)
    su_real = get_umbra(loc, t2)

    app.vars['d_est'] = '{:.7}'.format(dm)
    app.vars['d_real'] = '{:.7}'.format(reald)
    app.vars['su_est'] = '{:.4}'.format(su_est)
    app.vars['su_real'] = '{:.4}'.format(su_real)
    app.vars['lat'] = '{:.4}'.format(loc.latitude.value)
    app.vars['long'] = '{:.4}'.format(loc.longitude.value)

    sun_alt = '{:.4}'.format(sun_altaz.alt.value)
    percent_diff = '{:.2}'.format((dm - reald) / reald * 100.)

    print(loc, t2, t3)
    print(app.vars['c2time'], app.vars['c3time'])
    print(request.form['c2time'], request.form['c3time'])
    print(su_est, app.vars['su_est'])

    return render_template('view.html', location=app.vars['location'],
                           c2time=app.vars['c2time'], d_est=app.vars['d_est'], d_real=app.vars['d_real'],
                           su_real=app.vars['su_real'], lat=app.vars['lat'], su_est=app.vars['su_est'],
                           long=app.vars['long'], sun_alt=sun_alt, percent_diff=percent_diff)


# Called when you click clear button
@app.route('/clear')
def app_clear():
    clear_values()
    return redirect('/index')


def clear_values():
    for key in app.vars.keys():
        app.vars[key] = ''
    return

