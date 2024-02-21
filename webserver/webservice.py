import logging
import flask
import json

def mainpage():
    return flask.render_template( 'fom.html' )


# In reality, you'd probably do this with a template
def hello_world():
    return """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Hello</title>
</head>
<body>
    <p>Hello, world.</p>
</body>
</html>
"""

def x2( number ):
    number = float( number )
    resp = flask.make_response( str(2*number) )
    resp.headers['Content-Type'] = "text/plain; charset=utf=8"
    return resp

def multanddiv2( number ):
    number = float( number )
    resp = flask.make_response( { 'number': number,
                                  '×2': number*2,
                                  '÷2': number/2 } )
    return resp
    
def create_app():
    app = flask.Flask( __name__, instance_relative_config=True )
    app.logger.setLevel( logging.INFO )

    app.add_url_rule( "/", view_func=mainpage )
    app.add_url_rule( "/hello", view_func=hello_world )
    app.add_url_rule( "/mult2/<number>", view_func=x2 )
    app.add_url_rule( "/multdiv2/<number>", view_func=multanddiv2 )
    
    return app


