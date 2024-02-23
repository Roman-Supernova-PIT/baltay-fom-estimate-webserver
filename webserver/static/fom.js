import { rkWebUtil } from "./rkwebutil.js";

// Namespace

var fom = {};

// **********************************************************************
// **********************************************************************
// **********************************************************************

fom.Context = function()
{
    
    this.params = {'length1': { 'default': 1.9, 'desc': "Years of Survey 1" },
                   'length2': { 'default': 1.9, 'desc': "Years of Survey 2" },
                   'length3': { 'default': 1.9, 'desc': "Years of Survey 3" },
                   'sqdegim1': { 'default': 19.04, 'desc': "Square degrees of Survey 1" },
                   'sqdegim2': { 'default': 4.20, 'desc': "Square degrees of Survey 2" },
                   'sqdegim3': { 'default': 0.0, 'desc': "Square Degrees of Survey 3" },
                   'tfixim1': { 'default': 115.0, 'desc': "Imaging exposure time (s) of Survey 1" },
                   'tfixim2': { 'default': 450.0, 'desc': "Imaging exposure time (s) of Survey 2" },
                   'tfixim3': { 'default': 0, 'desc': "Imaging exposure time (s) of Survey 3" },
                   'z1': { 'default': 2.0, 'desc': "Redshift limit of Survey 1" },
                   'z2': { 'default': 2.0, 'desc': "Redshift limit of Survey 2" },
                   'z3': { 'default': 0.0, 'desc': "Redshift limit of Survey 3" },
                    // 'nlim': { 'default': 20, 'desc': "nlim(?)" },
                   'eff': { 'default': 0.9, 'desc': "Efficiency (frac of SNe caught)" },
                   'constsysim': { 'default': 0.015, 'desc': "Imaging systematic error at redshift divim" },
                   'divim': { 'default': 1.8, 'desc': "1+z where systematic error is constsysim" },
                   'stnim': { 'default': 10.0, 'desc': "S/N for imaging" },
                   'sqdegspec1': { 'default': 3.36, 'desc': "Square degrees for spectroscopy, Survey 1" },
                   'sqdegspec2': { 'default': 1.12, 'desc': "Square degrees for spectroscopy, Survey 2" },
                   'sqdegspec3': { 'default': 0.0, 'desc': "Square degrees for spectroscopy, Survey 3" },
                   'tfixspec1': { 'default': 900.0, 'desc': "Spectroscopy exposure time (s), Survey 1" },
                   'tfixspec2': { 'default': 3600.0, 'desc': "Spectroscopy exposure time (s), Survey 2" },
                   'tfixspec3': { 'default': 0, 'desc': "Spectroscopy exposure time (s), Survey 3" },
                   'constsysspec': { 'default': 0.01, 'desc': "Spectroscopy systematic error at redshit divspec" },
                   'divspec': { 'default': 1.8, 'desc': "1+z at which systematic error is constsysspec" },
                   'stnspec': { 'default': 20.0, 'desc': "Some sort of spectroscopic S/N..." },
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

