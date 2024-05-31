import { rkWebUtil } from "./rkwebutil.js";

// Namespace

var fom = {};

// **********************************************************************
// **********************************************************************
// **********************************************************************

fom.Context = function()
{
    
    this.params = {'length1': { 'default': 1.9, 'desc': "Years of Survey tier 1 (wide)" },
                   'length2': { 'default': 1.9, 'desc': "Years of Survey tier 2 (medium)" },
                   'length3': { 'default': 1.9, 'desc': "Years of Survey tier 3 (deep)" },
                   'sqdegim1': { 'default': 19.04, 'desc': "Square degrees of Survey tier 1 (wide)" },
                   'sqdegim2': { 'default': 4.20, 'desc': "Square degrees of Survey tier 2 (medium)" },
                   'sqdegim3': { 'default': 0.0, 'desc': "Square Degrees of Survey tier 3 (deep)" },
                   'tfixim1': { 'default': 115.0, 'desc': "Imaging exposure time (s) of Survey tier 1 (wide)" },
                   'tfixim2': { 'default': 450.0, 'desc': "Imaging exposure time (s) of Survey tier 2 (medium)" },
                   'tfixim3': { 'default': 0, 'desc': "Imaging exposure time (s) of Survey tier 3 (deep)" },
                   'z1': { 'default': 2.0, 'desc': "Imaging redshift limit of Survey tier 1 (wide)" },
                   'z2': { 'default': 2.0, 'desc': "Imaging redshift limit of Survey tier 2 (medium)" },
                   'z3': { 'default': 0.0, 'desc': "Imaging redshift limit of Survey tier 3 (deep)" },
                    // 'nlim': { 'default': 20, 'desc': "nlim(?)" },
                   'eff': { 'default': 0.9, 'desc': "Efficiency (frac of SNe caught)" },
                   'constsysim': { 'default': 0.015, 'desc': "Imaging systematic error coefficient" },
                   'divim': { 'default': 1.8, 'desc': "Imaging sysematic error denominator (1+z where systematic error is above)" },
                   'stnim': { 'default': 10.0, 'desc': "S/N for imaging to accept event" },
                   'sqdegspec1': { 'default': 3.36, 'desc': "Square degrees for spectroscopy, Survey tier 1 (wide)" },
                   'sqdegspec2': { 'default': 1.12, 'desc': "Square degrees for spectroscopy, Survey tier 2 (medium)" },
                   'sqdegspec3': { 'default': 0.0, 'desc': "Square degrees for spectroscopy, Survey tier 3 (deep)" },
                   'tfixspec1': { 'default': 900.0, 'desc': "Spectroscopy exposure time (s), Survey tier 1 (wide)" },
                   'tfixspec2': { 'default': 3600.0, 'desc': "Spectroscopy exposure time (s), Survey tier 2 (medium)" },
                   'tfixspec3': { 'default': 0, 'desc': "Spectroscopy exposure time (s), Survey tier 3 (deep)" },
                   'z1spec': { 'default': 2.0, 'desc': "Spectrscopy redshift limit of Survey tier 1 (wide)" },
                   'z2spec': { 'default': 2.0, 'desc': "Spectrscopy redshift limit of Survey tier 2 (medium)" },
                   'z3spec': { 'default': 0.0, 'desc': "Spectrscopy redshift limit of Survey tier 3 (deep)" },
                   'constsysspec': { 'default': 0.01, 'desc': "Spectroscopy systematic error coefficient" },
                   'divspec': { 'default': 1.8, 'desc': "spectroscopy sysematic error denominator (1+z where systematic error is above)" },
                   'stnspec': { 'default': 20.0, 'desc': "S/N for spectroscopy to accept event" },
                   'detail': { 'default': 0, 'desc': "Detailed output {0,1}" }
                  };
    this.paramwidgets = {};
};

fom.Context.prototype.renderpage = function()
{
    var self = this;
    
    this.maindiv = document.getElementById("pagebody");

    var table, tr, td;
    var table = rkWebUtil.elemaker( "table", this.maindiv );
    for ( let param in this.params ) {
        tr = rkWebUtil.elemaker( "tr", table );
        td = rkWebUtil.elemaker( "td", tr );
        this.paramwidgets[ param ] = rkWebUtil.elemaker( "input", td,
                                                         { "attributes": { "type": "text",
                                                                           "size": 8,
                                                                           "value": this.params[param].default
                                                                         } } );
        td = rkWebUtil.elemaker( "td", tr, { "text": this.params[param].desc } );
    }

    let button = rkWebUtil.elemaker( "button", this.maindiv, { "text": "Submit",
                                                               "click": function() { self.submit() } } )

    this.resultsdiv = rkWebUtil.elemaker( "div", this.maindiv );
};

fom.Context.prototype.submit = function()
{
    var self = this;

    rkWebUtil.wipeDiv( this.resultsdiv );

    var data = {};
    for ( let param in this.params ) {
        data[ param ] = parseFloat( this.paramwidgets[ param ].value );
    }
    
    let connector = new rkWebUtil.Connector( "/" );
    connector.sendHttpRequest( "calculate", data,
                               function( resp ) { self.showresults( resp ); } );
}

fom.Context.prototype.showresults = function( data )
{
    rkWebUtil.elemaker( "pre", this.resultsdiv, { "text": data.text } );
}
    
    


export { fom }

