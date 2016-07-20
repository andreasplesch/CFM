/** @namespace x3dom.nodeTypes */
/*
 * X3DOM JavaScript Library
 * http://www.x3dom.org
 *
 * (C)2009 Fraunhofer IGD, Darmstadt, Germany
 * Dual licensed under the MIT and GPL
 */

/* ### GeoCoordinate ### */
x3dom.registerNodeType(
    "GeoCoordinate",
    "Geospatial",
    defineClass(x3dom.nodeTypes.X3DCoordinateNode,
        
        /**
         * Constructor for GeoCoordinate
         * @constructs x3dom.nodeTypes.GeoCoordinate
         * @x3d 3.3
         * @component Geospatial
         * @status full
         * @extends x3dom.nodeTypes.X3DCoordinateNode
         * @param {Object} [ctx=null] - context object, containing initial settings like namespace
         * @classdesc The GeoCoordinate node specifies a list of coordinates in a spatial reference frame.
         * It is used in the coord field of vertex-based geometry nodes including IndexedFaceSet, IndexedLineSet, and PointSet.
         */
        function (ctx) {
            x3dom.nodeTypes.GeoCoordinate.superClass.call(this, ctx);


            /**
             * The point array is used to contain the actual geospatial coordinates and should be provided in a format consistent with that specified for the particular geoSystem.
             * @var {x3dom.fields.MFVec3f} point
             * @memberof x3dom.nodeTypes.GeoCoordinate
             * @initvalue []
             * @field x3dom
             * @instance
             */
            this.addField_MFVec3f(ctx, 'point', []);

            /**
             * The geoSystem field is used to define the spatial reference frame.
             * @var {x3dom.fields.MFString} geoSystem
             * @range {["GD", ...], ["UTM", ...], ["GC", ...]}
             * @memberof x3dom.nodeTypes.GeoCoordinate
             * @initvalue ['GD','WE']
             * @field x3dom
             * @instance
             */
            this.addField_MFString(ctx, 'geoSystem', ['GD', 'WE']);

            /**
             * The geoOrigin field is used to specify a local coordinate frame for extended precision.
             * @var {x3dom.fields.SFNode} geoOrigin
             * @memberof x3dom.nodeTypes.GeoCoordinate
             * @initvalue x3dom.nodeTypes.GeoOrigin
             * @field x3dom
             * @instance
             */
            this.addField_SFNode('geoOrigin', x3dom.nodeTypes.GeoOrigin);
        
        },
        {
            elipsoideParameters:
            {
                'AA' : [ 'Airy 1830', '6377563.396', '299.3249646' ],
                'AM' : [ 'Modified Airy', '6377340.189', '299.3249646' ],
                'AN' : [ 'Australian National', '6378160', '298.25' ],
                'BN' : [ 'Bessel 1841 (Namibia)', '6377483.865', '299.1528128' ],
                'BR' : [ 'Bessel 1841 (Ethiopia Indonesia...)', '6377397.155', '299.1528128' ],
                'CC' : [ 'Clarke 1866', '6378206.4', '294.9786982' ],
                'CD' : [ 'Clarke 1880', '6378249.145', '293.465' ],
                'EA' : [ 'Everest (India 1830)', '6377276.345', '300.8017' ],
                'EB' : [ 'Everest (Sabah & Sarawak)', '6377298.556', '300.8017' ],
                'EC' : [ 'Everest (India 1956)', '6377301.243', '300.8017' ],
                'ED' : [ 'Everest (W. Malaysia 1969)', '6377295.664', '300.8017' ],
                'EE' : [ 'Everest (W. Malaysia & Singapore 1948)', '6377304.063', '300.8017' ],
                'EF' : [ 'Everest (Pakistan)', '6377309.613', '300.8017' ],
                'FA' : [ 'Modified Fischer 1960', '6378155', '298.3' ],
                'HE' : [ 'Helmert 1906', '6378200', '298.3' ],
                'HO' : [ 'Hough 1960', '6378270', '297' ],
                'ID' : [ 'Indonesian 1974', '6378160', '298.247' ],
                'IN' : [ 'International 1924', '6378388', '297' ],
                'KA' : [ 'Krassovsky 1940', '6378245', '298.3' ],
                'RF' : [ 'Geodetic Reference System 1980 (GRS 80)', '6378137', '298.257222101' ],
                'SA' : [ 'South American 1969', '6378160', '298.25' ],
                'WD' : [ 'WGS 72', '6378135', '298.26' ],
                'WE' : [ 'WGS 84', '6378137', '298.257223563' ]
            },

            fieldChanged: function(fieldName) {
                if (fieldName == "point" || fieldName == "geoSystem") {
                    Array.forEach(this._parentNodes, function (node) {
                        node.fieldChanged("coord");
                    });
                }
            },

            isLogitudeFirst: function(geoSystem) {
                for(var i=0; i<geoSystem.length; ++i)
                    if(geoSystem[i] == 'longitude_first')
                        return true;

                return false;
            },

            getElipsoideCode: function(geoSystem)
            {
                for(var i=0; i<geoSystem.length; ++i)
                {
                    var code = geoSystem[i];
                    if(this.elipsoideParameters[code])
                        return code;
                }
                //default elipsoide code
                return 'WE';
            },

            getElipsoide: function(geoSystem)
            {
                return this.elipsoideParameters[this.getElipsoideCode(geoSystem)];
            },

            getReferenceFrame: function(geoSystem)
            {
                for(var i=0; i<geoSystem.length; ++i)
                {
                    var code = geoSystem[i];

                    if(code == 'GD' || code == 'GDC')
                        return 'GD';
                    if(code == 'GC' || code == 'GCC')
                        return 'GC';
                    if(code == 'UTM')
                        return 'UTM';

                    else
                        x3dom.debug.logError('Unknown GEO system: [' + geoSystem + ']');
                }

                // default reference frame is GD WE
                return 'GD';
            },

            getUTMZone: function(geoSystem)
            {
                for(var i=0; i<geoSystem.length; ++i)
                {
                    var code = geoSystem[i];

                    if(code[0] == 'Z')
                        return code.substring(1);
                }
                // no zone found
                x3dom.debug.logError('no UTM zone but is required:' + geoSystem);
            },

            getUTMHemisphere: function(geoSystem)
            {
                for(var i=0; i<geoSystem.length; ++i)
                {
                    var code = geoSystem[i];

                    if(code == 'S')
                        return code;
                }
                // North by default according to spec
                return 'N';
            },

            isUTMEastingFirst: function(geoSystem)
            {
                for(var i=0; i<geoSystem.length; ++i)
                {
                    var code = geoSystem[i];
                    if(code == 'easting_first')
                        return true;
                }
                // Northing first by default according to spec
                return false;
            },

            UTMtoGC: function(geoSystem, coords)
            {
                //parse UTM projection parameters
                var utmzone = this.getUTMZone(geoSystem);
                if(utmzone < 1 || utmzone > 60 || utmzone === undefined)
                    return x3dom.debug.logError('invalid UTM zone: ' + utmzone + ' in geosystem ' + geoSystem);
                var hemisphere = this.getUTMHemisphere(geoSystem);
                var eastingFirst = this.isUTMEastingFirst(geoSystem);
                var elipsoide = this.getElipsoide(geoSystem);
                //below from U.W. Green Bay Prof. Dutch; returns coordinates in the input ell., not WGS84
                var a = elipsoide[1];
                var f = 1/elipsoide[2];
                var k0 = 0.9996; //scale on central meridian
                var b = a * (1 - f); //polar axis.
                var esq = (1 - (b/a)*(b/a)); //e squared for use in expansions
                var e = Math.sqrt(esq); //eccentricity
                var e0 = e/Math.sqrt(1 - esq); //Called e prime in reference
                var e0sq = esq/(1 - esq); // e0 squared - always even powers
                var zcm = 3 + 6 * (utmzone - 1) - 180; //Central meridian of zone
                var e1 = (1 - Math.sqrt(1 - esq))/(1 + Math.sqrt(1 - esq)); //Called e1 in USGS PP 1395 also
                var e1sq = e1*e1;
                //var M0 = 0; //In case origin other than zero lat - not needed for standard UTM
                var output = new x3dom.fields.MFVec3f();
                var rad2deg = 180/Math.PI;

                var f3o64 = 3/64;
                var f5o256 = 5/256;
                var f27o32 = 27/32;
                var f21o16 = 21/16;
                var f55o32 = 55/32;
                var f151o96 = 151/96;
                var f1097o512 = 1097/512;
                var fmua = 1 - esq*(0.25 + esq*(f3o64 + f5o256*esq));
                var fphi11 = e1*(1.5 - f27o32*e1sq);
                var fphi12 = e1sq*(f21o16 - f55o32*e1sq);
                var current, x, y, z, M, mu, phi1, cosphi1, C1, tanphi1, T1, T1sq;
                var esinphi1, oneesinphi1, N1, R1, D, Dsq, C1sq, phi, lng;
                
                for(var i=0; i<coords.length; ++i)
                {
                    x = (eastingFirst ? coords[i].x : coords[i].y);
                    y = (eastingFirst ? coords[i].y : coords[i].x);
                    z = coords[i].z;

                    //var M = M0 + y/k0; //Arc length along standard meridian.
                    //var M = y/k0;
                    //if (hemisphere == "S"){ M = M0 + (y - 10000000)/k; }
                    M = (hemisphere == "S" ? (y - 10000000) : y )/k0 ;
                    mu = M/(a * fmua);
                    phi1 = mu + fphi11*Math.sin(2*mu) + fphi12*Math.sin(4*mu); //Footprint Latitude
                    phi1 = phi1 + e1*(e1sq*(Math.sin(6*mu)*f151o96 + Math.sin(8*mu)*f1097o512));
                    cosphi1 = Math.cos(phi1);
                    C1 = e0sq*cosphi1*cosphi1;
                    tanphi1 = Math.tan(phi1);
                    T1 = tanphi1*tanphi1;
                    T1sq = T1*T1;
                    esinphi1 = e*Math.sin(phi1);
                    oneesinphi1 = 1 - esinphi1*esinphi1;
                    N1 = a/Math.sqrt(oneesinphi1);
                    //R1 = N1*(1-e*e)/oneesinphi1;
                    R1 = N1*(1-esq)/oneesinphi1;
                    D = (x-500000)/(N1*k0);
                    Dsq = D*D;
                    C1sq = C1*C1;
                    phi = Dsq*(0.5 - Dsq*(5 + 3*T1 + 10*C1 - 4*C1sq - 9*e0sq)/24);
                    phi = phi + Math.pow(D,6)*(61 + 90*T1 + 298*C1 + 45*T1sq -252*e0sq - 3*C1sq)/720;
                    phi = phi1 - (N1*tanphi1/R1)*phi;
                    lng = D*(1 + Dsq*((-1 -2*T1 -C1)/6 + Dsq*(5 - 2*C1 + 28*T1 - 3*C1sq +8*e0sq + 24*T1sq)/120))/cosphi1;
                    current = new x3dom.fields.SFVec3f();
                    current.x = zcm + rad2deg*lng;
                    current.y = rad2deg*phi;
                    current.z = coords[i].z;
                    output.push(current);
                }
                //x3dom.debug.logInfo('transformed coords ' + output);

                //GD to GC and return
                var GDgeoSystem = new x3dom.fields.MFString();
                // there may be a better way to construct this geoSystem
                GDgeoSystem.push("GD");
                GDgeoSystem.push(this.getElipsoideCode(geoSystem));
                GDgeoSystem.push("longitude_first");
                return this.GDtoGC(GDgeoSystem, output);
            },
            
            //just used for GeoPositionInterpolator now; after slerp
            //for coordinates in the same UTM zone stripe, linear interpolation is almost identical
            //so this is for correctness and if UTM zone is used outside stripe
            GCtoUTM: function(geoSystem, coords) {
                // geoSystem is target UTM
                // GCtoGD
                var coordsGD = this.GCtoGD(geoSystem, coords);
                // GDtoUTM
                // parse UTM projection parameters
                var utmzone = this.getUTMZone(geoSystem);
                if(utmzone < 1 || utmzone > 60 || utmzone === undefined)
                    return x3dom.debug.logError('invalid UTM zone: ' + utmzone + ' in geosystem ' + geoSystem);
                var hemisphere = this.getUTMHemisphere(geoSystem);
                var eastingFirst = this.isUTMEastingFirst(geoSystem);
                var elipsoide = this.getElipsoide(geoSystem);
                //below from U.W. Green Bay Prof. Dutch; returns coordinates in the input ell., not WGS84
                var a = elipsoide[1];
                var f = 1/elipsoide[2];
                var k0 = 0.9996; //scale on central meridian
                var b = a * (1 - f); //polar axis.
                var esq = (1 - (b/a)*(b/a)); //e squared for use in expansions
                var e = Math.sqrt(esq); //eccentricity
                //var e0 = e/Math.sqrt(1 - esq); //Called e prime in reference
                var e0sq = esq/(1 - esq); // e0 squared - always even powers
                //var e1 = (1 - Math.sqrt(1 - esq))/(1 + Math.sqrt(1 - esq)); //Called e1 in USGS PP 1395 also
                //var e1sq = e1*e1;
                var M0 = 0; //In case origin other than zero lat - not needed for standard UTM
                var deg2rad = Math.PI/180;
                var zcmrad = (3 + 6 * (utmzone - 1) - 180)*deg2rad; //Central meridian of zone
                var coordsUTM = new x3dom.fields.MFVec3f();
                var N, T, C, A, M, x, y, phi, lng, cosphi, tanphi, Asq;
                var i, current;
                var fMphi = 1 - esq*(1/4 + esq*(3/64 + 5*esq/256));
                var fM2phi = esq*(3/8 + esq*(3/32 + 45*esq/1024));
                var fM4phi = esq*esq*(15/256 + esq*45/1024);
                var fM6phi = esq*esq*esq*(35/3072);
                for (i=0; i<coordsGD.length; ++i) {
                    current = new x3dom.fields.SFVec3f();
                    phi = coordsGD[i].y*deg2rad;
                    lng = coordsGD[i].x*deg2rad;
                    cosphi = Math.cos(phi);
                    tanphi = Math.tan(phi);
                    
                    N = a/Math.sqrt(1 - Math.pow(e * Math.sin(phi), 2));
                    T = Math.pow(tanphi, 2);
                    C = e0sq*Math.pow(cosphi, 2);
                    A = (lng - zcmrad) * cosphi;
                    //Calculate M
                    M = phi*fMphi;
                    M = M - Math.sin(2*phi)*fM2phi;
                    M = M + Math.sin(4*phi)*fM4phi;
                    M = M - Math.sin(6*phi)*fM6phi;
                    M = M * a;//Arc length along standard meridian
                    //Calculate UTM Values
                    Asq = A*A;
                    x = k0*N*A*(1 + Asq*((1 - T + C)/6 + Asq*(5 - T*(18 + T) + 72*C - 58*e0sq)/120));//Easting relative to CM
                    x = x + 500000;//Easting standard 
                    y = k0*(M - M0 + N*tanphi*(Asq*(0.5 + Asq*((5 - T + 9*C + 4*C*C)/24 + Asq*(61 - T*(58 + T) + 600*C - 330*e0sq)/720))));//Northing from equator
                    if (y < 0) {
                        if (hemisphere == "N") {
                            x3dom.debug.logError('UTM zone in northern hemisphere but coordinates in southern!');
                        }
                        y = 10000000+y;}
                    current.x = eastingFirst ? x : y;
                    current.y = eastingFirst ? y : x;
                    current.z = coordsGD[i].z;
                    coordsUTM.push(current);
                }
                return coordsUTM;    
            },

            GDtoGC: function(geoSystem, coords) {

                var output = new x3dom.fields.MFVec3f();

                var elipsoide = this.getElipsoide(geoSystem);
                //var radius = elipsoide[1];
                var A = elipsoide[1];
                var eccentricity = elipsoide[2];

                var longitudeFirst = this.isLogitudeFirst(geoSystem);

                // large parts of this code from freeWRL
                var A2 = A*A;
                var F = 1.0/eccentricity;
                var C = A*(1.0-F);
                //var C2 = C*C;
                var C2oA2 = C*C/A2;
                var Eps2 = F*(2.0-F);
                //var Eps25 = 0.25*Eps2;

                var radiansPerDegree = 0.0174532925199432957692;

                // for (current in coords)
                for(var i=0; i<coords.length; ++i)
                {
                    var current = new x3dom.fields.SFVec3f();

                    var source_lat = radiansPerDegree * (longitudeFirst == true ? coords[i].y : coords[i].x);
                    var source_lon = radiansPerDegree * (longitudeFirst == true ? coords[i].x : coords[i].y);

                    var slat = Math.sin(source_lat);
                    var slat2 = slat*slat;
                    var clat = Math.cos(source_lat);

                    /* square root approximation for Rn */
                    /* replaced by real sqrt
                     var Rn = A / ( (0.25 - Eps25 * slat2 + 0.9999944354799/4.0) +
                     (0.25-Eps25 * slat2)/(0.25 - Eps25 * slat2 + 0.9999944354799/4.0));
                     */

                    // with real sqrt; really slower ?
                    var Rn = A / Math.sqrt(1.0 - Eps2 * slat2);

                    var RnPh = Rn + coords[i].z;

                    current.x = RnPh * clat * Math.cos(source_lon);
                    current.y = RnPh * clat * Math.sin(source_lon);
                    current.z = (C2oA2 * Rn + coords[i].z) * slat;

                    output.push(current);
                }

                return output;
            },
            
            GCtoGD: function(geoSystem, coords) {
                // below uses http://info.ogp.org.uk/geodesy/guides/docs/G7-2.pdf
                var output = new x3dom.fields.MFVec3f();
                var rad2deg = 180/Math.PI;
                var ellipsoide = this.getElipsoide(geoSystem);
                var a = ellipsoide[1];
                var f = 1/ellipsoide[2];
                var b = a * (1 - f); //polar axis.
                var esq = (1 - (b/a)*(b/a)); //e squared for use in expansions
                //var e = Math.sqrt(esq); //eccentricity
                var eps = esq/(1 - esq);
                var x, y, z, p, q, lat, nu, elev, lon;
                var current;
                for(var i=0; i<coords.length; ++i) {
                    current = new x3dom.fields.SFVec3f();
                    x = coords[i].x;
                    y = coords[i].y;
                    z = coords[i].z;
                    
                    p = Math.sqrt(x*x + y*y);
                    q = Math.atan((z * a) / (p * b));
                    lat = Math.atan((z + eps * b * Math.pow(Math.sin(q),3))/(p - esq * a * Math.pow(Math.cos(q),3)));
                    nu = a / Math.sqrt(1-esq * Math.pow(Math.sin(lat),2));
                    elev = p/Math.cos(lat) - nu;
                    // atan2 gets the sign correct for longitude; is exact since in circular section
                    lon = Math.atan2(y, x);
                    
                    current.x = lon * rad2deg;
                    current.y = lat * rad2deg;
                    current.z = elev;

                    output.push(current);
                }
                return output;
            },
        
            GEOtoGC: function(geoSystem, geoOrigin, coords) {
                
                var referenceFrame = this.getReferenceFrame(geoSystem);

                if(referenceFrame == 'GD') {
                    return this.GDtoGC(geoSystem, coords);
                }

                else if(referenceFrame == 'UTM') {
                    return this.UTMtoGC(geoSystem, coords);
                }

                else if(referenceFrame ==  'GC')
                {
                    // Performance Hack
                    // Normaly GDtoGC & UTMtoGC will create a copy
                    // If we are already in GC & have an origin: we have to copy here
                    // Else Origin will change original DOM elements

                    if(geoOrigin.node)
                    {
                        var copy = new x3dom.fields.MFVec3f();
                        for(var i=0; i<coords.length; ++i)
                        {
                            var current = new x3dom.fields.SFVec3f();

                            current.x = coords[i].x;
                            current.y = coords[i].y;
                            current.z = coords[i].z;

                            copy.push(current);
                        }
                        return copy;
                    }
                    else
                        return coords;
                }
                else {
                    x3dom.debug.logError('Unknown geoSystem: ' + geoSystem[0]);
                    return new x3dom.fields.MFVec3f();
                }
            },

            GCtoGEO: function(geoSystem, geoOrigin, coords) {
                
                var referenceFrame = this.getReferenceFrame(geoSystem);

                if(referenceFrame == 'GD') {
                    var coordsGD = this.GCtoGD(geoSystem, coords);
                    if(!this.isLogitudeFirst(geoSystem)) {
                        var currentx;
                        for(var i=0; i<coordsGD.length; ++i) {
                            currentx = coordsGD[i].x;
                            coordsGD[i].x = coordsGD[i].y;
                            coordsGD[i].y = currentx;
                        }
                    }
                    return coordsGD;
                }

                else if(referenceFrame == 'UTM')
                    return this.GCtoUTM(geoSystem, coords);

                else if(referenceFrame ==  'GC')
                {
                    // Performance Hack
                    // Normaly GDtoGC & UTMtoGC will create a copy
                    // If we are already in GC & have an origin: we have to copy here
                    // Else Origin will change original DOM elements

                    if(geoOrigin.node)
                    {
                        var copy = new x3dom.fields.MFVec3f();
                        for(var i=0; i<coords.length; ++i)
                        {
                            var current = new x3dom.fields.SFVec3f();

                            current = coords[i].copy;
                            //current.y = coords[i].y;
                            //current.z = coords[i].z;

                            copy.push(current);
                        }
                        return copy;
                    }
                    else
                        return coords;
                }
                else {
                    x3dom.debug.logError('Unknown geoSystem: ' + geoSystem[0]);
                    return new x3dom.fields.MFVec3f();
                }
            },
            

            OriginToGC: function(geoOrigin)
            {
                // dummy function to send a scalar to an array function
                var geoCoords = geoOrigin.node._vf.geoCoords;
                var geoSystem = geoOrigin.node._vf.geoSystem;

                var point = new x3dom.fields.SFVec3f();
                point.x = geoCoords.x;
                point.y = geoCoords.y;
                point.z = geoCoords.z;

                var temp = new x3dom.fields.MFVec3f();
                temp.push(point);

                // transform origin to GeoCentric
                var origin = this.GEOtoGC(geoSystem, geoOrigin, temp);

                return origin[0];
            },
            
            GCtoX3D: function(geoSystem, geoOrigin, coords)
            {
                // needs geoSystem since GC can have different ellipsoids
                var gc = coords;

                // transform by origin
                if(geoOrigin.node)
                {
                    // transform points by origin
                    var origin = this.OriginToGC(geoOrigin);
                    //avoid expensive matrix inversion by just subtracting origin
                    var matrix = x3dom.fields.SFMatrix4f.translation(origin.negate());
                    
                    //also rotate Y up if requested
                    if(geoOrigin.node._vf.rotateYUp)
                    {
                        //rotation is inverse of GeoLocation rotation, eg. Up to Y and N to -Z
                        var rotmat = x3dom.nodeTypes.GeoLocation.prototype.getGeoRotMat(geoSystem, origin).inverse();
                        //first translate, then rotate
                        matrix = rotmat.mult(matrix);
                    }
                    
                    for(var i=0; i<coords.length; ++i)
                        gc[i] = matrix.multMatrixPnt(gc[i]);
                }

                return gc;
            },

            GEOtoX3D: function(geoSystem, geoOrigin, coords)
            {
                // transform points to GeoCentric
                var gc = this.GEOtoGC(geoSystem, geoOrigin, coords);
                return this.GCtoX3D(geoSystem, geoOrigin, gc);
                
                // transform by origin
                if(geoOrigin.node)
                {
                    // transform points by origin
                    var origin = this.OriginToGC(geoOrigin);
                    //avoid expensive matrix inversion by just subtracting origin
                    var matrix = x3dom.fields.SFMatrix4f.translation(origin.negate());
                    
                    //also rotate Y up if requested
                    if(geoOrigin.node._vf.rotateYUp)
                    {
                        //rotation is inverse of GeoLocation rotation, eg. Up to Y and N to -Z
                        var rotmat = x3dom.nodeTypes.GeoLocation.prototype.getGeoRotMat(geoSystem, origin).inverse();
                        //first translate, then rotate
                        matrix = rotmat.mult(matrix);
                    }
                    
                    for(var i=0; i<coords.length; ++i)
                        gc[i] = matrix.multMatrixPnt(gc[i]);
                }

                return gc;
            },

            getPoints: function()
            {
                return this.GEOtoX3D(this._vf.geoSystem, this._cf.geoOrigin, this._vf.point);
            }
        }
    )
);
/** @namespace x3dom.nodeTypes */
/*
 * X3DOM JavaScript Library
 * http://www.x3dom.org
 *
 * (C)2009 Fraunhofer IGD, Darmstadt, Germany
 * Dual licensed under the MIT and GPL
 */

/* ### GeoElevationGrid ### */
x3dom.registerNodeType(
    "GeoElevationGrid",
    "Geospatial",
    defineClass(x3dom.nodeTypes.X3DGeometryNode,
        
        /**
         * Constructor for GeoElevationGrid
         * @constructs x3dom.nodeTypes.GeoElevationGrid
         * @x3d 3.3
         * @component Geospatial
         * @status experimental
         * @extends x3dom.nodeTypes.X3DGeometryNode
         * @param {Object} [ctx=null] - context object, containing initial settings like namespace
         * @classdesc The GeoElevationGrid node specifies a uniform grid of elevation values within some spatial reference frame.
         * These are then transparently transformed into a geocentric, curved-earth representation.
         */
        function (ctx) {
            x3dom.nodeTypes.GeoElevationGrid.superClass.call(this, ctx);


            /**
             * The geoSystem field is used to define the spatial reference frame.
             * @var {x3dom.fields.MFString} geoSystem
             * @range {["GD", ...], ["UTM", ...], ["GC", ...]}
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue ['GD','WE']
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'geoSystem', ['GD', 'WE']);

            /**
             * The geoGridOrigin field specifies the geographic coordinate for the south-west corner (bottom-left) of the dataset.
             * @var {x3dom.fields.SFVec3d} geoGridOrigin
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue 0,0,0
             * @field x3d
             * @instance
             */
            this.addField_SFVec3d(ctx, 'geoGridOrigin', 0, 0, 0);

            /**
             * The height array contains xDimension × zDimension floating point values that represent elevation above the ellipsoid or the geoid, as appropriate.
             * These values are given in row-major order from west to east, south to north.
             * @var {x3dom.fields.MFDouble} height
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue 0,0
             * @field x3d
             * @instance
             */
            this.addField_MFDouble(ctx, 'height', 0, 0);

            /**
             * The ccw field defines the ordering of the vertex coordinates of the geometry with respect to user-given or automatically generated normal vectors used in the lighting model equations.
             * @var {x3dom.fields.SFBool} ccw
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue true
             * @field x3d
             * @instance
             */
            this.addField_SFBool(ctx, 'ccw', true);
            //this.addField_SFBool(ctx, 'colorPerVertex', true);

            /**
             * The creaseAngle field affects how default normals are generated.
             * If the angle between the geometric normals of two adjacent faces is less than the crease angle, normals shall be calculated so that the faces are shaded smoothly across the edge; otherwise, normals shall be calculated so that a lighting discontinuity across the edge is produced.
             * Crease angles shall be greater than or equal to 0.0 angle base units.
             * @var {x3dom.fields.SFDouble} creaseAngle
             * @range [0, inf]
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue 0
             * @field x3d
             * @instance
             */
            this.addField_SFDouble(ctx, 'creaseAngle', 0);
            //this.addField_SFBool(ctx, 'normalPerVertex', true);
            //this.addField_SFBool(ctx, 'solid', true);

            /**
             * Defines the grid size in x.
             * @var {x3dom.fields.SFInt32} xDimension
             * @range [0, inf]
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue 0
             * @field x3d
             * @instance
             */
            this.addField_SFInt32(ctx, 'xDimension', 0);

            /**
             * When the geoSystem is "GD", xSpacing refers to the number of units of longitude in angle base units between adjacent height values.
             * When the geoSystem is "UTM", xSpacing refers to the number of eastings (length base units) between adjacent height values
             * @var {x3dom.fields.SFDouble} xSpacing
             * @range [0, inf]
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue 1.0
             * @field x3d
             * @instance
             */
            this.addField_SFDouble(ctx, 'xSpacing', 1.0);

            /**
             * The yScale value can be used to produce a vertical exaggeration of the data when it is displayed. If this value is set greater than 1.0, all heights will appear larger than actual.
             * @var {x3dom.fields.SFFloat} yScale
             * @range [0, inf]
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue 1
             * @field x3d
             * @instance
             */
            this.addField_SFFloat(ctx, 'yScale', 1);

            /**
             * Defines the grid size in z.
             * @var {x3dom.fields.SFInt32} zDimension
             * @range [0, inf]
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue 0
             * @field x3d
             * @instance
             */
            this.addField_SFInt32(ctx, 'zDimension', 0);

            /**
             * When the geoSystem is "GD", zSpacing refers to the number of units of latitude in angle base units between vertical height values.
             * When the geoSystem is "UTM", zSpacing refers to the number of northings (length base units) between vertical height values.
             * @var {x3dom.fields.SFDouble} zSpacing
             * @range [0, inf]
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue 1.0
             * @field x3d
             * @instance
             */
            this.addField_SFDouble(ctx, 'zSpacing', 1.0);
            // this.addField_SFNode('color', x3dom.nodeTypes.PropertySetGeometry);
            // this.addField_SFNode('normal', x3dom.nodeTypes.PropertySetGeometry);
            // this.addField_SFNode('texCoord', x3dom.nodeTypes.PropertySetGeometry);

            /**
             * The geoOrigin field is used to specify a local coordinate frame for extended precision.
             * @var {x3dom.fields.SFNode} geoOrigin
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue x3dom.nodeTypes.GeoOrigin
             * @field x3d
             * @instance
             */
            this.addField_SFNode('geoOrigin', x3dom.nodeTypes.GeoOrigin);

            /**
             * Specifies whether this geometry should be rendered with or without lighting.
             * @var {x3dom.fields.SFBool} lit
             * @memberof x3dom.nodeTypes.GeoElevationGrid
             * @initvalue true
             * @field x3dom
             * @instance
             */
            this.addField_SFBool(ctx, 'lit', true);
        
        },
        {
            nodeChanged: function()
            {
                var geoSystem = this._vf.geoSystem;
                var geoOrigin = this._cf.geoOrigin;

                var height = this._vf.height;

                var yScale = this._vf.yScale;
                var xDimension = this._vf.xDimension;
                var zDimension = this._vf.zDimension;
                var xSpacing = this._vf.xSpacing;
                var zSpacing = this._vf.zSpacing;
                var geoGridOrigin = this._vf.geoGridOrigin;

                // check for no height == dimensions
                if(height.length !== (xDimension * zDimension))
                    x3dom.debug.logError('GeoElevationGrid: height.length(' + height.length +
                        ') != x/zDimension(' + xDimension + '*' + zDimension + ')');

                var longitude_first = x3dom.nodeTypes.GeoCoordinate.prototype.isLogitudeFirst(geoSystem);
                var ccw = this._vf.ccw;

                // coords, texture coords
                var delta_x = 1 / (xDimension-1);
                var delta_z = 1 / (zDimension-1);

                var positions = new x3dom.fields.MFVec3f();
                var texCoords = new x3dom.fields.MFVec2f();

                for(var z=0; z<zDimension; ++z)
                    for(var x=0; x<xDimension; ++x)
                    {
                        // texture coord
                        var tex_coord = new x3dom.fields.SFVec2f(x*delta_x, z*delta_z);
                        texCoords.push(tex_coord);

                        // coord
                        var coord = new x3dom.fields.SFVec3f();
                        if(longitude_first)
                        {
                            coord.x = x * xSpacing;
                            coord.y = z * zSpacing;
                        }
                        else
                        {
                            coord.x = z * zSpacing;
                            coord.y = x * xSpacing;
                        }
                        coord.z = height[(z*xDimension)+x] * yScale;
                        coord = coord.add(geoGridOrigin);

                        positions.push(coord);
                    }

                // indices
                var indices = new x3dom.fields.MFInt32();
                for(var z=0; z<(zDimension-1); z++)
                {
                    for(var x=0; x<(xDimension-1); x++)
                    {
                        var p0 = x + (z * xDimension);
                        var p1 = x + (z * xDimension) + 1;
                        var p2 = x + ((z + 1) * xDimension) + 1;
                        var p3 = x + ((z + 1) * xDimension);

                        if(ccw)
                        {
                            indices.push(p0);
                            indices.push(p1);
                            indices.push(p2);

                            indices.push(p0);
                            indices.push(p2);
                            indices.push(p3);
                        }
                        else
                        {
                            indices.push(p0);
                            indices.push(p3);
                            indices.push(p2);

                            indices.push(p0);
                            indices.push(p2);
                            indices.push(p1);
                        }
                    }
                }

                // convert to x3dom coord system
                var transformed = x3dom.nodeTypes.GeoCoordinate.prototype.GEOtoX3D(geoSystem, geoOrigin, positions);

                //if we want flat shading, we have to duplicate some vertices here
                //(as webgl does only support single-indexed rendering)
                if (this._vf.creaseAngle <= x3dom.fields.Eps) {

                    var that = this;

                    (function (){
                        var indicesFlat   = new x3dom.fields.MFInt32(),
                            positionsFlat = new x3dom.fields.MFVec3f(),
                            texCoordsFlat = new x3dom.fields.MFVec3f();

                        that.generateNonIndexedTriangleData(indices, transformed, null, texCoords, null,
                            positionsFlat, null, texCoordsFlat, null);

                        for (var i = 0; i < positionsFlat.length; ++i) {
                            indicesFlat.push(i);
                        }

                        that._mesh._indices[0]   = indicesFlat.toGL();
                        that._mesh._positions[0] = positionsFlat.toGL();
                        that._mesh._texCoords[0] = texCoordsFlat.toGL();
                    })();

                    this._mesh.calcNormals(0);
                }
                //smooth shading
                else {
                    this._mesh._indices[0]   = indices.toGL();
                    this._mesh._positions[0] = transformed.toGL();
                    this._mesh._texCoords[0] = texCoords.toGL();

                    this._mesh.calcNormals(Math.PI);
                }

                this._mesh._invalidate = true;
                this._mesh._numFaces = this._mesh._indices[0].length / 3;
                this._mesh._numCoords = this._mesh._positions[0].length / 3;
            },

            generateNonIndexedTriangleData: function(indices, positions, normals, texCoords, colors,
                                                     newPositions, newNormals, newTexCoords, newColors)
            {
                //@todo: add support for RGBA colors and 3D texture coordinates
                //@todo: if there is any need for that, add multi-index support

                for (var i = 0; i < indices.length; i+=3) {
                    var i0 = indices[i  ],
                        i1 = indices[i+1],
                        i2 = indices[i+2];

                    if (positions) {
                        var p0 = new x3dom.fields.SFVec3f(),
                            p1 = new x3dom.fields.SFVec3f(),
                            p2 = new x3dom.fields.SFVec3f();

                        p0.setValues(positions[i0]);
                        p1.setValues(positions[i1]);
                        p2.setValues(positions[i2]);

                        newPositions.push(p0);
                        newPositions.push(p1);
                        newPositions.push(p2);
                    }

                    if (normals) {
                        var n0 = new x3dom.fields.SFVec3f(),
                            n1 = new x3dom.fields.SFVec3f(),
                            n2 = new x3dom.fields.SFVec3f();

                        n0.setValues(normals[i0]);
                        n1.setValues(normals[i1]);
                        n2.setValues(normals[i2]);

                        newNormals.push(n0);
                        newNormals.push(n1);
                        newNormals.push(n2);
                    }

                    if (texCoords) {
                        var t0 = new x3dom.fields.SFVec2f(),
                            t1 = new x3dom.fields.SFVec2f(),
                            t2 = new x3dom.fields.SFVec2f();

                        t0.setValues(texCoords[i0]);
                        t1.setValues(texCoords[i1]);
                        t1.setValues(texCoords[i2]);

                        newTexCoords.push(t0);
                        newTexCoords.push(t1);
                        newTexCoords.push(t2);
                    }

                    if (colors) {
                        var c0 = new x3dom.fields.SFVec3f(),
                            c1 = new x3dom.fields.SFVec3f(),
                            c2 = new x3dom.fields.SFVec3f();

                        c0.setValues(texCoords[i0]);
                        c1.setValues(texCoords[i1]);
                        c1.setValues(texCoords[i2]);

                        newColors.push(c0);
                        newColors.push(c1);
                        newColors.push(c2);
                    }
                }
            }
        }
    )
);/** @namespace x3dom.nodeTypes */
/*
 * X3DOM JavaScript Library
 * http://www.x3dom.org
 *
 * (C)2009 Fraunhofer IGD, Darmstadt, Germany
 * Dual licensed under the MIT and GPL
 */

/* ### GeoLOD ### */
x3dom.registerNodeType(
    "GeoLOD",
    "Geospatial",
    defineClass(x3dom.nodeTypes.X3DLODNode,
        
        /**
         * Constructor for GeoLOD
         * @constructs x3dom.nodeTypes.GeoLOD
         * @x3d 3.3
         * @component Geospatial
         * @status experimental
         * @extends x3dom.nodeTypes.X3DLODNode
         * @param {Object} [ctx=null] - context object, containing initial settings like namespace
         * @classdesc The GeoLOD node provides a terrain-specialized form of the LOD node.
         * It is a grouping node that specifies two different levels of detail for an object using a tree structure, where 0 to 4 children can be specified, and also efficiently manages the loading and unloading of these levels of detail.
         */
        function (ctx) {
            x3dom.nodeTypes.GeoLOD.superClass.call(this, ctx);


            /**
             * The geoSystem field is used to define the spatial reference frame.
             * @var {x3dom.fields.MFString} geoSystem
             * @range {["GD", ...], ["UTM", ...], ["GC", ...]}
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue ['GD','WE']
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'geoSystem', ['GD', 'WE']);

            /**
             * The rootUrl and rootNode fields provide two different ways to specify the geometry of the root tile.
             * You may use one or the other. The rootNode field lets you include the geometry for the root tile directly within the X3D file.
             * @var {x3dom.fields.MFString} rootUrl
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue []
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'rootUrl', []);

            /**
             * When the viewer enters the specified range, this geometry is replaced with the contents of the four children files defined by child1Url through child4Url.
             * The four children files are loaded into memory only when the user is within the specified range. Similarly, these are unloaded from memory when the user leaves this range.
             * @var {x3dom.fields.MFString} child1Url
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue []
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'child1Url', []);

            /**
             * When the viewer enters the specified range, this geometry is replaced with the contents of the four children files defined by child1Url through child4Url.
             * The four children files are loaded into memory only when the user is within the specified range. Similarly, these are unloaded from memory when the user leaves this range.
             * @var {x3dom.fields.MFString} child2Url
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue []
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'child2Url', []);

            /**
             * When the viewer enters the specified range, this geometry is replaced with the contents of the four children files defined by child1Url through child4Url.
             * The four children files are loaded into memory only when the user is within the specified range. Similarly, these are unloaded from memory when the user leaves this range.
             * @var {x3dom.fields.MFString} child3Url
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue []
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'child3Url', []);

            /**
             * When the viewer enters the specified range, this geometry is replaced with the contents of the four children files defined by child1Url through child4Url.
             * The four children files are loaded into memory only when the user is within the specified range. Similarly, these are unloaded from memory when the user leaves this range.
             * @var {x3dom.fields.MFString} child4Url
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue []
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'child4Url', []);
            //this.addField_SFVec3d(ctx, 'center', 0, 0, 0);

            /**
             * The level of detail is switched depending upon whether the user is closer or farther than range length base units from the geospatial coordinate center.
             * @var {x3dom.fields.SFFloat} range
             * @range [0, inf]
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue 10
             * @field x3d
             * @instance
             */
            this.addField_SFFloat(ctx, 'range', 10);

            /**
             *
             * @var {x3dom.fields.SFString} referenceBindableDescription
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue []
             * @field x3dom
             * @instance
             */
            this.addField_SFString(ctx, 'referenceBindableDescription', []);

            /**
             * The geoOrigin field is used to specify a local coordinate frame for extended precision.
             * @var {x3dom.fields.SFNode} geoOrigin
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue x3dom.nodeTypes.X3DChildNode
             * @field x3d
             * @instance
             */
            this.addField_SFNode('geoOrigin', x3dom.nodeTypes.X3DChildNode);

            /**
             * The rootUrl and rootNode fields provide two different ways to specify the geometry of the root tile. The rootUrl field lets you specify a URL for a file that contains the geometry.
             * @var {x3dom.fields.SFNode} rootNode
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue x3dom.nodeTypes.X3DChildNode
             * @field x3dom
             * @instance
             */
            this.addField_SFNode('rootNode', x3dom.nodeTypes.X3DChildNode);

            /**
             *
             * @var {x3dom.fields.SFNode} privateChild1Node
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue x3dom.nodeTypes.X3DChildNode
             * @field x3dom
             * @instance
             */
            this.addField_SFNode('privateChild1Node', x3dom.nodeTypes.X3DChildNode);

            /**
             *
             * @var {x3dom.fields.SFNode} privateChild2Node
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue x3dom.nodeTypes.X3DChildNode
             * @field x3dom
             * @instance
             */
            this.addField_SFNode('privateChild2Node', x3dom.nodeTypes.X3DChildNode);

            /**
             *
             * @var {x3dom.fields.SFNode} privateChild3Node
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue x3dom.nodeTypes.X3DChildNode
             * @field x3dom
             * @instance
             */
            this.addField_SFNode('privateChild3Node', x3dom.nodeTypes.X3DChildNode);

            /**
             *
             * @var {x3dom.fields.SFNode} privateChild4Node
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue x3dom.nodeTypes.X3DChildNode
             * @field x3dom
             * @instance
             */
            this.addField_SFNode('privateChild4Node', x3dom.nodeTypes.X3DChildNode);

            /**
             *
             * @var {x3dom.fields.SFNode} privateRootNode
             * @memberof x3dom.nodeTypes.GeoLOD
             * @initvalue x3dom.nodeTypes.X3DChildNode
             * @field x3dom
             * @instance
             */
            this.addField_SFNode('privateRootNode', x3dom.nodeTypes.X3DChildNode);
        
        }
    )
);/** @namespace x3dom.nodeTypes */
/*
 * X3DOM JavaScript Library
 * http://www.x3dom.org
 *
 * (C)2009 Fraunhofer IGD, Darmstadt, Germany
 * Dual licensed under the MIT and GPL
 */

/* ### GeoLocation ### */
x3dom.registerNodeType(
    "GeoLocation",
    "Geospatial",
    //was X3DGroupingNode which is how the node is defined in the spec
    defineClass(x3dom.nodeTypes.X3DTransformNode,
        
        /**
         * Constructor for GeoLocation
         * @constructs x3dom.nodeTypes.GeoLocation
         * @x3d 3.3
         * @component Geospatial
         * @status experimental
         * @extends x3dom.nodeTypes.X3DTransformNode
         * @param {Object} [ctx=null] - context object, containing initial settings like namespace
         * @classdesc The GeoLocation node provides the ability to geo-reference any standard X3D model.
         * That is, to take an ordinary X3D model, contained within the children of the node, and to specify its geospatial location.
         * This node is a grouping node that can be thought of as a Transform node.
         * However, the GeoLocation node specifies an absolute location, not a relative one, so content developers should not nest GeoLocation nodes within each other.
         */
        function (ctx) {
            x3dom.nodeTypes.GeoLocation.superClass.call(this, ctx);


            /**
             * The geoSystem field is used to define the spatial reference frame.
             * @var {x3dom.fields.MFString} geoSystem
             * @range {["GD", ...], ["UTM", ...], ["GC", ...]}
             * @memberof x3dom.nodeTypes.GeoLocation
             * @initvalue ['GD','WE']
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'geoSystem', ['GD', 'WE']);

            /**
             * The geometry of the nodes in children is to be specified in units of metres in X3D coordinates relative to the location specified by the geoCoords field.
             * The geoCoords field can be used to dynamically update the geospatial location of the model.
             * @var {x3dom.fields.SFVec3d} geoCoords
             * @memberof x3dom.nodeTypes.GeoLocation
             * @initvalue 0,0,0
             * @field x3d
             * @instance
             */
            this.addField_SFVec3d(ctx, 'geoCoords', 0, 0, 0);

            /**
             * The geoOrigin field is used to specify a local coordinate frame for extended precision.
             * @var {x3dom.fields.SFNode} geoOrigin
             * @memberof x3dom.nodeTypes.GeoLocation
             * @initvalue x3dom.nodeTypes.GeoOrigin
             * @field x3d
             * @instance
             */
            this.addField_SFNode('geoOrigin', x3dom.nodeTypes.GeoOrigin);
        },

        {
            nodeChanged: function()
            {
                // similar to what transform in Grouping.js does
                var position = this._vf.geoCoords;
                var geoSystem = this._vf.geoSystem;
                var geoOrigin = this._cf.geoOrigin; // gets only populated if in nodeChanged()

                this._trafo =  this.getGeoTransRotMat(geoSystem, geoOrigin, position);
            },
        
            getGeoRotMat: function (geoSystem, positionGC)
            {
                //returns transformation matrix to align coordinate system with geoposition as required:
                //2 rotations to get required orientation
                //Up (Y) to skywards, and depth (-Z) to North
                //1) around X to point up by
                //angle between Z and new up plus 90
                //(angle between Z and orig. up)
                //2) around Z to get orig. up on longitude

                //var newUp = positionGC.normalize();
                var coords = new x3dom.fields.MFVec3f();
                coords.push(positionGC);
                var positionGD = x3dom.nodeTypes.GeoCoordinate.prototype.GCtoGD(geoSystem, coords)[0];
                
                var Xaxis = new  x3dom.fields.SFVec3f(1,0,0);

                // below uses geocentric latitude but only geodetic latitude would give exact tangential plane
                // http://info.ogp.org.uk/geodesy/guides/docs/G7-2.pdf
                // has formulas for deriving geodetic latitude, eg a GCtoGD function
                //var rotlat = Math.PI - Math.asin(newUp.z); // latitude as asin of z; only valid for spheres
                var rotlat = 180 - positionGD.y; // latitude
                var deg2rad = Math.PI/180;
                var rotUpQuat = x3dom.fields.Quaternion.axisAngle(Xaxis, rotlat*deg2rad);

                //var rotlon = Math.PI/2 + Math.atan2(newUp.y, newUp.x);// 90 to get to prime meridian; atan2 gets the sign correct for longitude; is exact since in circular section
                var rotlon = 90 + positionGD.x;// 90 to get to prime meridian; atan2 gets the sign correct for longitude; is exact since in circular section
                var Zaxis = new x3dom.fields.SFVec3f(0,0,1);
                var rotZQuat = x3dom.fields.Quaternion.axisAngle(Zaxis, rotlon*deg2rad);

                //return rotZQuat.toMatrix().mult(rotUpQuat.toMatrix();
                return rotZQuat.multiply(rotUpQuat).toMatrix();

            },

            getGeoTransRotMat: function (geoSystem, geoOrigin, position)
            {
                // accept geocoords, return translation/rotation transform matrix
                var coords = new x3dom.fields.MFVec3f();
                coords.push(position);

                var transformed = x3dom.nodeTypes.GeoCoordinate.prototype.GEOtoGC(geoSystem, geoOrigin, coords)[0];
                var rotMat = this.getGeoRotMat(geoSystem, transformed);

                // account for geoOrigin with and without rotateYUp
                if (geoOrigin.node)
                {
                    var origin = x3dom.nodeTypes.GeoCoordinate.prototype.OriginToGC(geoOrigin);
                    if(geoOrigin.node._vf.rotateYUp)
                    {
                        // inverse rotation after original rotation and offset
                        // just skipping all rotations produces incorrect position
                        var rotMatOrigin = this.getGeoRotMat(geoSystem, origin);
                        return rotMatOrigin.inverse().mult(x3dom.fields.SFMatrix4f.translation(transformed.subtract(origin)).mult(rotMat));
                    }
                    //rotate, then translate; account for geoOrigin by subtracting origin from GeoLocation
                    return x3dom.fields.SFMatrix4f.translation(transformed.subtract(origin)).mult(rotMat);
                }
                else
                //no GeoOrigin: first rotate, then translate
                {
                    return x3dom.fields.SFMatrix4f.translation(transformed).mult(rotMat);
                }
            },

            //mimic what transform node does
            fieldChanged: function (fieldName)
            {
                if (fieldName == "geoSystem" || fieldName == "geoCoords" ||
                    fieldName == "geoOrigin")
                {
                    var position = this._vf.geoCoords;
                    var geoSystem = this._vf.geoSystem;
                    var geoOrigin = this._cf.geoOrigin;
                    this._trafo =  this.getGeoTransRotMat(geoSystem, geoOrigin, position);

                    this.invalidateVolume();
                    //this.invalidateCache();
                }
                else if (fieldName == "render") {
                    this.invalidateVolume();
                    //this.invalidateCache();
                }
            }
            //deal with geolocation in geolocation here? behaviour is undefined in spec

        }
    )
);
/** @namespace x3dom.nodeTypes */
/*
 * X3DOM JavaScript Library
 * http://www.x3dom.org
 *
 * (C)2009 Fraunhofer IGD, Darmstadt, Germany
 * Dual licensed under the MIT and GPL
 */

/* ### GeoMetadata ### */
x3dom.registerNodeType(
    "GeoMetadata",
    "Geospatial",
    defineClass(x3dom.nodeTypes.X3DInfoNode,
        
        /**
         * Constructor for GeoMetadata
         * @constructs x3dom.nodeTypes.GeoMetadata
         * @x3d 3.3
         * @component Geospatial
         * @status full
         * @extends x3dom.nodeTypes.X3DInfoNode
         * @param {Object} [ctx=null] - context object, containing initial settings like namespace
         * @classdesc The GeoMetadata node supports the specification of metadata describing any number of geospatial nodes.
         * This is similar to a WorldInfo node, but specifically for describing geospatial information.
         */
        function (ctx) {
            x3dom.nodeTypes.GeoMetadata.superClass.call(this, ctx);


            /**
             * The url field is used to specify a hypertext link to an external, complete metadata description.
             * Multiple URL strings can be specified in order to provide alternative locations for the same metadata information.
             * The summary field may be used to specify the format of the metadata in the case where this cannot be deduced easily.
             * @var {x3dom.fields.MFString} url
             * @memberof x3dom.nodeTypes.GeoMetadata
             * @initvalue []
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'url', []);

            /**
             * The data field is used to list all of the other nodes in a scene by DEF name that reference the data described in the GeoMetadata node.
             * The nodes in the data field are not rendered, so DEF/USE can be used in order to first describe them and then to use them in the scene graph.
             * @var {x3dom.fields.MFNode} data
             * @memberof x3dom.nodeTypes.GeoMetadata
             * @initvalue x3dom.nodeTypes.X3DInfoNode
             * @field x3d
             * @instance
             */
            this.addField_MFNode('data', x3dom.nodeTypes.X3DInfoNode);

            /**
             * The summary string array contains a set of keyword/value pairs, with each keyword and its subsequent value contained in a separate string; i.e., there should always be an even (or zero) number of strings.
             * @var {x3dom.fields.MFString} summary
             * @memberof x3dom.nodeTypes.GeoMetadata
             * @initvalue []
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'summary', []);
        
        }
    )
);/** @namespace x3dom.nodeTypes */
/*
 * X3DOM JavaScript Library
 * http://www.x3dom.org
 *
 * (C)2009 Fraunhofer IGD, Darmstadt, Germany
 * Dual licensed under the MIT and GPL
 */

/* ### GeoOrigin ### */
x3dom.registerNodeType(
    "GeoOrigin",
    "Geospatial",
    defineClass(x3dom.nodeTypes.X3DNode,
        
        /**
         * Constructor for GeoOrigin
         * @constructs x3dom.nodeTypes.GeoOrigin
         * @x3d 3.2
         * @component Geospatial
         * @status full
         * @extends x3dom.nodeTypes.X3DNode
         * @param {Object} [ctx=null] - context object, containing initial settings like namespace
         * @classdesc The GeoOrigin node defines an absolute geospatial location and an implicit local coordinate frame against which geometry is referenced.
         * This node is used to translate from geographical coordinates into a local Cartesian coordinate system which can be managed by the X3D browser. This node is deprecated as of X3D 3.3
         */
        function (ctx) {
            x3dom.nodeTypes.GeoOrigin.superClass.call(this, ctx);


            /**
             * The geoSystem field is used to define the spatial reference frame.
             * @var {x3dom.fields.MFString} geoSystem
             * @range {["GD", ...], ["UTM", ...], ["GC", ...]}
             * @memberof x3dom.nodeTypes.GeoOrigin
             * @initvalue ['GD','WE']
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'geoSystem', ['GD', 'WE']);

            /**
             * The geoCoords field is used to specify a local coordinate frame for extended precision.
             * @var {x3dom.fields.SFVec3d} geoCoords
             * @memberof x3dom.nodeTypes.GeoOrigin
             * @initvalue 0,0,0
             * @field x3d
             * @instance
             */
            this.addField_SFVec3d(ctx, 'geoCoords', 0, 0, 0);

            /**
             * The rotateYUp field is used to specify whether coordinates of nodes that use this GeoOrigin are to be rotated such that their up direction is aligned with the X3D Y axis.
             * The default behavior is to not perform this operation.
             * @var {x3dom.fields.SFBool} rotateYUp
             * @memberof x3dom.nodeTypes.GeoOrigin
             * @initvalue false
             * @field x3d
             * @instance
             */
            this.addField_SFBool(ctx, 'rotateYUp', false);
        
        }
    )
);/** @namespace x3dom.nodeTypes */
/*
 * X3DOM JavaScript Library
 * http://www.x3dom.org
 *
 * (C)2009 Fraunhofer IGD, Darmstadt, Germany
 * Dual licensed under the MIT and GPL
 */

/* ### GeoPositionInterpolator ### */
x3dom.registerNodeType(
    "GeoPositionInterpolator",
    "Geospatial",
    defineClass(x3dom.nodeTypes.X3DInterpolatorNode,
        
        /**
         * Constructor for GeoPositionInterpolator
         * @constructs x3dom.nodeTypes.GeoPositionInterpolator
         * @x3d 3.3
         * @component Geospatial
         * @status experimental
         * @extends x3dom.nodeTypes.X3DInterpolatorNode
         * @param {Object} [ctx=null] - context object, containing initial settings like namespace
         * @classdesc The GeoPositionInterpolator node provides an interpolator capability where key values are specified in geographic coordinates and the interpolation is performed within the specified spatial reference frame.
         */
        function (ctx) {
            x3dom.nodeTypes.GeoPositionInterpolator.superClass.call(this, ctx);


            /**
             * The geoSystem field is used to define the spatial reference frame.
             * @var {x3dom.fields.MFString} geoSystem
             * @range {["GD", ...], ["UTM", ...], ["GC", ...]}
             * @memberof x3dom.nodeTypes.GeoPositionInterpolator
             * @initvalue ['GD','WE']
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'geoSystem', ['GD', 'WE']);

            /**
             * The keyValue array is used to contain the actual coordinates and should be provided in a format consistent with that specified for the particular geoSystem.
             * @var {x3dom.fields.MFVec3d} keyValue
             * @memberof x3dom.nodeTypes.GeoPositionInterpolator
             * @initvalue []
             * @field x3d
             * @instance
             */
            this.addField_MFVec3f(ctx, 'keyValue', []);

            /**
             * The geoOrigin field is used to specify a local coordinate frame for extended precision.
             * @var {x3dom.fields.SFNode} geoOrigin
             * @memberof x3dom.nodeTypes.GeoPositionInterpolator
             * @initvalue x3dom.nodeTypes.X3DInterpolatorNode
             * @field x3d
             * @instance
             */
            this.addField_SFNode('geoOrigin', x3dom.nodeTypes.X3DInterpolatorNode);
            
            /**
             * The onGreatCircle field is used to specify whether coordinates will be interpolated along a great circle path
             * The default behavior is to not perform this operation.
             * @var {x3dom.fields.SFBool} onGreatCircle
             * @memberof x3dom.nodeTypes.GeoPositionInterpolator
             * @initvalue false
             * @field x3dom
             * @instance
             */
            this.addField_SFBool(ctx, 'onGreatCircle', false);
            
            /**
             * The exactGreatCircle field is used to specify wether the exact Great Circle path is used or a linear approximation which avoids coordinate transformations
             * The default behavior is to perform this operation.
             * @var {x3dom.fields.SFBool} onGreatCircle
             * @memberof x3dom.nodeTypes.GeoPositionInterpolator
             * @initvalue true
             * @field x3dom
             * @instance
             */
            this.addField_SFBool(ctx, 'exactGreatCircle', false);
        
            /*
             * original simply does linear interpolation, then converts to geocentric for value_changed
             */
            
            /* optional: onGreatCircle=true
             * use slerp directly on geocentric vectors. slerp is in fields.js for quaternions, just set 4th comp. to zero
             * should work ok, elevation should be also interpolated ok with the slerp of the gc vectors
             * then just convert to geosystem (implement GCtoUTM)
             */
             
        },
        {
            // adapted from X3DInterpolator.js
            
            // we need our own key and keyValue, uses hint for larger arrays and returns also found index
            linearInterpHintKeyValue: function (time, keyHint, key, keyValue, interp) {
                // add guess as input where to find time; should be close to index of last time?
                // do wraparound search in both directions since interpolation often goes back
                var keylength = key.length;
                
                if (time <= key[0])
                    return [0, keyValue[0]];
                
                else if (time >= key[keylength - 1])
                    return [keylength - 1, keyValue[keylength - 1]];
                
                var keyIndexStart = keyHint ;
                var i;
                // strictly loop only to keylength/2 but does not hurt
                for (var ii = 0; ii < keylength - 1; ++ii) {
                    // look forward
                    i = (keyIndexStart + ii) % keylength;
                    //i+1 can be outside array but undefined leads to false in check
                    if ((key[i] < time) && (time <= key[i+1]))
                        return [i, interp( keyValue[i], keyValue[i+1],
                                (time - key[i]) / (key[i+1] - key[i]) )];
                    // look backward
                    i = (keyIndexStart - ii + keylength) % keylength; 
                    if ((key[i] < time) && (time <= key[i+1]))
                        return [i, interp( keyValue[i], keyValue[i+1],
                                (time - key[i]) / (key[i+1] - key[i]) )];                    
                }
                return [0, keyValue[0]];
            },
            
            // adapted from fields.js
            slerp: function (a, b, t) {
                // calculate the cosine
                // since a and b are not unit vectors here; this is the only real change
                var cosom = a.dot(b)/(a.length()*b.length());
                var rot1;
            
                /* 
                 * does not apply for geometric slerp
                 * adjust signs if necessary
                if (cosom < 0.0)
                {
                    cosom = -cosom;
                    rot1 = b.negate();
                }
                else
                */
                {
                    rot1 = new x3dom.fields.SFVec3f(b.x, b.y, b.z);
                }
            
                // calculate interpolating coeffs
                var scalerot0, scalerot1;
                
                if ((1.0 - cosom) > 0.00001)
                {
                    // standard case
                    var omega = Math.acos(cosom);
                    var sinom = Math.sin(omega);
                    scalerot0 = Math.sin((1.0 - t) * omega) / sinom;
                    scalerot1 = Math.sin(t * omega) / sinom;
                }
                else
                {
                    // rot0 and rot1 very close - just do linear interp.
                    scalerot0 = 1.0 - t;
                    scalerot1 = t;
                }
            
                // build the new vector
                return a.multiply(scalerot0).add(rot1.multiply(scalerot1));
            },
            
            nodeChanged: function() {
                // set up initial values
                this._keyValueGC = x3dom.nodeTypes.GeoCoordinate.prototype.GEOtoGC(this._vf.geoSystem, this._vf.geoOrigin, this._vf.keyValue);
                this._keyHint = 0;
                // sanity check key.length vs. keyValue.length
                
                // "linearize" great circle path to minimize coordinate conversions
                // an angular segment size of 0.1 degrees produces a misfit of 2.5m vertically, eg. good enough
                // eg. max. arc is 180, 1800 segments, seems to much
                // lets say max. 180 segments per interval, so for 18 degrees arc still good fit
                // larger 18 degrees just use a 180 equal segments: 180/180=1 degrees max, prod. ca. 250m misfit 
                // go through each interval
                // measure angular size
                // if size < 18 degrees
                //   divide interval equally such segments are close to 0.1:
                //     number of new segments is int(size/0.1)+1 or so; 0.45/0.1=4.5=4, plus 1;ok
                // otherwise number of new segments is 180
                // produce new keys, scale by original key
                // produce corresponding new key values: GC on great circle
                var a, b, t, cosom, n_segments, omega, newKey, keyMin, keyMax, keyDist;
                var maxangle = 18 * Math.PI/180 ;
                var maxsegment = 0.1 * Math.PI/180 ;
                //var cosmax_angle = Math.cos(maxangle);
                var n_segments_max = maxangle/maxsegment;
                this._linKey = new x3dom.fields.MFFloat();
                for (var i = 0; i < this._vf.key.length-1; ++i) {
                    a = this._keyValueGC[i];
                    b = this._keyValueGC[i+1];
                    cosom = a.dot(b)/(a.length()*b.length());
                    omega = Math.acos(cosom);
                    n_segments = omega < maxangle ? Math.floor(omega/maxsegment)+1 : n_segments_max;
                    keyMin = this._vf.key[i];
                    keyMax = this._vf.key[i+1];
                    keyDist = keyMax - keyMin;
                    for (var n = 0; n < n_segments; ++n) {
                        newKey = keyMin + (n/n_segments)*keyDist;
                        this._linKey.push(newKey);
                    }                    
                }
                // still need last key
                this._linKey.push(keyMax);
                // go through all new keys as times and produce new keyValues in given geoSystem
                this._linKeyValue = new x3dom.fields.MFVec3f();
                var indexLinValueGC, linvalue, coords;
                var hint = 0;
                for (i = 0; i < this._linKey.length; ++i) {
                    indexLinValueGC = this.linearInterpHintKeyValue(this._linKey[i], hint, this._vf.key, this._keyValueGC, x3dom.nodeTypes.GeoPositionInterpolator.prototype.slerp);
                    hint = indexLinValueGC[0];
                    coords = new x3dom.fields.MFVec3f();
                    coords.push(indexLinValueGC[1]);
                    linvalue = x3dom.nodeTypes.GeoCoordinate.prototype.GCtoGEO(this._vf.geoSystem, this._vf.geoOrigin, coords)[0];
                    this._linKeyValue.push(linvalue);
                }
                //need to treat GD specially if going across date line                    
                var referenceFrame = x3dom.nodeTypes.GeoCoordinate.prototype.getReferenceFrame(this._vf.geoSystem);
                if(referenceFrame == 'GD') {
                    var isLongitudeFirst = x3dom.nodeTypes.GeoCoordinate.prototype.isLogitudeFirst(this._vf.geoSystem);
                    var val1, valMid, val3, lon1, lonMid, lon3, nlonMid, extrap;
                    //does not work if going across in very last segment ...
                    for (i = 0; i < this._linKey.length - 2; ++i) {
                            val1 = this._linKeyValue[i] ;                        
                            valMid = this._linKeyValue[i+1] ;
                            val3 = this._linKeyValue[i+2] ;
                            lon1 = isLongitudeFirst ? val1.x : val1.y; 
                            lonMid = isLongitudeFirst ? valMid.x : valMid.y; 
                            lon3 = isLongitudeFirst ? val3.x : val3.y;
                            if (Math.abs(lon3-lon1) > 180) {
                                x3dom.debug.logError("date line: " + lon1 + " " + lonMid + " " + lon3);
                                //push lonMid to 180, and lon3 also but on the other side
                                nlonMid=lon1 > 0 ? 180 : -180;
                                //adjust Key position
                                extrap = (nlonMid-lon1)/(lonMid-lon1);
                                this._linKey[i+1] = this._linKey[i] + extrap*(this._linKey[i+1]-this._linKey[i]);
                                this._linKey[i+2] = this._linKey[i+1] ; // + 0.000000000000001;
                                //extrapolate xyz location, 
                                this._linKeyValue[i+1] = val1.add((valMid.subtract(val1)).multiply(extrap));
                                this._linKeyValue[i+2] = this._linKeyValue[i+1].copy();
                                //third point on other side
                                if (isLongitudeFirst) {
                                    this._linKeyValue[i+2].x = -nlonMid;
                                }
                                else {
                                    this._linKeyValue[i+2].y = -nlonMid;
                                }
                                //skip next
                                i = i + 1;              
                            }
                    }
                }
            },
            
            // adapted from PositionInterpolator.js
            fieldChanged: function(fieldName)
            {
                if(fieldName === "set_fraction")
                {
                    var value, indexValue, valueGC, valueX3D, coords ;
                    if(this._vf.onGreatCircle) {
                        if(this._vf.exactGreatCircle) {
                            indexValue = this.linearInterpHintKeyValue(this._vf.set_fraction, this._keyHint, this._vf.key, this._keyValueGC, x3dom.nodeTypes.GeoPositionInterpolator.prototype.slerp);
                            this._keyHint = indexValue[0];
                            valueGC = indexValue[1];                            
                            coords = new x3dom.fields.MFVec3f();
                            coords.push(valueGC);
                            value = x3dom.nodeTypes.GeoCoordinate.prototype.GCtoGEO(this._vf.geoSystem, this._vf.geoOrigin, coords)[0];
                        }
                        else {
                            // use linearized Great Circle
                            indexValue = this.linearInterpHintKeyValue(this._vf.set_fraction, this._keyHint, this._linKey, this._linKeyValue, function (a, b, t) {
                            return a.multiply(1.0-t).add(b.multiply(t));
                            });
                            this._keyHint = indexValue[0];
                            value = indexValue[1];
                            // just do slerp on orginal for GC
                            // do not reuse this._keyHint since it is a different index
                            // should be a small array
                            valueGC = this.linearInterpHintKeyValue(this._vf.set_fraction, 0, this._vf.key, this._keyValueGC, x3dom.nodeTypes.GeoPositionInterpolator.prototype.slerp);
                        }
                    }
                    else {
                        value = this.linearInterp(this._vf.set_fraction, function (a, b, t) {
                            return a.multiply(1.0-t).add(b.multiply(t));                        
                        });
                        coords = new x3dom.fields.MFVec3f();
                        coords.push(value);
                        valueGC = x3dom.nodeTypes.GeoCoordinate.prototype.GEOtoGC(this._vf.geoSystem, this._vf.geoOrigin, coords)[0];
                    }
                    //x3dom.debug.logInfo("interpolated fraction: " + this._vf.set_fraction);
                    
                    //x3dom.debug.logInfo("interpolated GD: " + value);
                    //x3dom.debug.logInfo("interpolated GC: " + this._keyValueGC);
                    //account for GeoOrigin, eg. transform to X3D coordinates
                    //coords.push(valueGC);
                    //valueX3D = x3dom.nodeTypes.GeoCoordinate.prototype.GCtoX3D(this._vf.geoSystem, this._vf.geoOrigin, coords)[0];
                    //this.postMessage('value_changed', valueX3D);
                    
                    this.postMessage('value_changed', valueGC);
                    this.postMessage('geovalue_changed', value);
                }
            }        
        
        }
    )
);

/** @namespace x3dom.nodeTypes */
/*
 * X3DOM JavaScript Library
 * http://www.x3dom.org
 *
 * (C)2009 Fraunhofer IGD, Darmstadt, Germany
 * Dual licensed under the MIT and GPL
 */

/* ### GeoTransform ### */
x3dom.registerNodeType(
    "GeoTransform",
    "Geospatial",
    defineClass(x3dom.nodeTypes.X3DGroupingNode,
        
        /**
         * Constructor for GeoTransform
         * @constructs x3dom.nodeTypes.GeoTransform
         * @x3d 3.3
         * @component Geospatial
         * @status full
         * @extends x3dom.nodeTypes.X3DGroupingNode
         * @param {Object} [ctx=null] - context object, containing initial settings like namespace
         * @classdesc The GeoTransform node is a grouping node that defines a coordinate system for its children to support the translation and orientation of geometry built using GeoCoordinate nodes within the local world coordinate system.
         * The X-Z plane of a GeoTransform coordinate system is tangent to the ellipsoid of the spatial reference frame at the location specified by the geoCenter field.
         */
        function (ctx) {
            x3dom.nodeTypes.GeoTransform.superClass.call(this, ctx);


            /**
             * The geoCenter field specifies, in the spatial reference frame specified by the geoSystem field, the location at which the local coordinate system is centered.
             * @var {x3dom.fields.SFVec3d} geoCenter
             * @memberof x3dom.nodeTypes.GeoTransform
             * @initvalue 0,0,0
             * @field x3d
             * @instance
             */
            this.addField_SFVec3d(ctx, 'geoCenter', 0, 0, 0);

            /**
             * Defines the rotation component of the transformation.
             * @var {x3dom.fields.SFRotation} rotation
             * @memberof x3dom.nodeTypes.GeoTransform
             * @initvalue 0,0,1,0
             * @field x3d
             * @instance
             */
            this.addField_SFRotation(ctx, 'rotation', 0, 0, 1, 0);

            /**
             * Defines the scale component of the transformation.
             * @var {x3dom.fields.SFVec3f} scale
             * @memberof x3dom.nodeTypes.GeoTransform
             * @initvalue 1,1,1
             * @field x3d
             * @instance
             */
            this.addField_SFVec3f(ctx, 'scale', 1, 1, 1);

            /**
             * The scaleOrientation specifies a rotation of the coordinate system before the scale (to specify scales in arbitrary orientations).
             * The scaleOrientation applies only to the scale operation.
             * @var {x3dom.fields.SFRotation} scaleOrientation
             * @memberof x3dom.nodeTypes.GeoTransform
             * @initvalue 0,0,1,0
             * @field x3d
             * @instance
             */
            this.addField_SFRotation(ctx, 'scaleOrientation', 0, 0, 1, 0);

            /**
             * The translation field specifies a translation to the coordinate system.
             * @var {x3dom.fields.SFVec3f} translation
             * @memberof x3dom.nodeTypes.GeoTransform
             * @initvalue 0,0,0
             * @field x3d
             * @instance
             */
            this.addField_SFVec3f(ctx, 'translation', 0, 0, 0);

            /**
             * The geoOrigin field is used to specify a local coordinate frame for extended precision.
             * @var {x3dom.fields.SFNode} geoOrigin
             * @memberof x3dom.nodeTypes.GeoTransform
             * @initvalue x3dom.nodeTypes.Transform
             * @field x3d
             * @instance
             */
            this.addField_SFNode('geoOrigin', x3dom.nodeTypes.Transform);

            /**
             * The geoSystem field is used to define the spatial reference frame.
             * @var {x3dom.fields.MFString} geoSystem
             * @range {["GD", ...], ["UTM", ...], ["GC", ...]}
             * @memberof x3dom.nodeTypes.GeoTransform
             * @initvalue ['GD','WE']
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'geoSystem', ['GD', 'WE']);
        
        }
    )
);
/** @namespace x3dom.nodeTypes */
/*
 * X3DOM JavaScript Library
 * http://www.x3dom.org
 *
 * (C)2009 Fraunhofer IGD, Darmstadt, Germany
 * Dual licensed under the MIT and GPL
 */

/* ### GeoViewpoint ### */
x3dom.registerNodeType(
    "GeoViewpoint",
    "Geospatial",
    defineClass(x3dom.nodeTypes.X3DViewpointNode,
        
        /**
         * Constructor for GeoViewpoint
         * @constructs x3dom.nodeTypes.GeoViewpoint
         * @x3d 3.3
         * @component Geospatial
         * @status experimental
         * @extends x3dom.nodeTypes.X3DViewpointNode
         * @param {Object} [ctx=null] - context object, containing initial settings like namespace
         * @classdesc The GeoViewpoint node allows the specification of a viewpoint in terms of a geospatial coordinate.
         * This node can be used wherever a Viewpoint node can be used and can be combined with Viewpoint nodes in the same scene.
         */
        function (ctx) {
            x3dom.nodeTypes.GeoViewpoint.superClass.call(this, ctx);


            /**
             * The geoSystem field is used to define the spatial reference frame.
             * @var {x3dom.fields.MFString} geoSystem
             * @range {["GD", ...], ["UTM", ...], ["GC", ...]}
             * @memberof x3dom.nodeTypes.GeoViewpoint
             * @initvalue ['GD','WE']
             * @field x3d
             * @instance
             */
            this.addField_MFString(ctx, 'geoSystem', ['GD', 'WE']);

            /**
             * Preferred minimum viewing angle from this viewpoint in radians.
             * Small field of view roughly corresponds to a telephoto lens, large field of view roughly corresponds to a wide-angle lens.
             * Hint: modifying Viewpoint distance to object may be better for zooming.
             * Warning: fieldOfView may not be correct for different window sizes and aspect ratios.
             * Interchange profile hint: this field may be ignored.
             * @var {x3dom.fields.SFFloat} fieldOfView
             * @range [0, pi]
             * @memberof x3dom.nodeTypes.GeoViewpoint
             * @initvalue 0.785398
             * @field x3d
             * @instance
             */
            this.addField_SFFloat(ctx, 'fieldOfView', 0.785398);

            /**
             * The orientation fields of the Viewpoint node specifies relative orientation to the default orientation.
             * @var {x3dom.fields.SFRotation} orientation
             * @memberof x3dom.nodeTypes.GeoViewpoint
             * @initvalue 0,0,1,0
             * @field x3d
             * @instance
             */
            this.addField_SFRotation(ctx, 'orientation', 0, 0, 1, 0);

            /**
             * The centerOfRotation field specifies a center about which to rotate the user's eyepoint when in EXAMINE mode.
             * The coordinates are provided in the coordinate system specified by geoSystem.
             * @var {x3dom.fields.SFVec3f} centerOfRotation
             * @memberof x3dom.nodeTypes.GeoViewpoint
             * @initvalue 0,0,0
             * @field x3d
             * @instance
             */
            this.addField_SFVec3f(ctx, 'centerOfRotation', 0, 0, 0);

            /**
             * The position fields of the Viewpoint node specifies a relative location in the local coordinate system.
             * The coordinates are provided in the coordinate system specified by geoSystem. 
             * @var {x3dom.fields.SFVec3d} position
             * @memberof x3dom.nodeTypes.GeoViewpoint
             * @initvalue 0,0,100000
             * @field x3d
             * @instance
             */
            this.addField_SFVec3d(ctx, 'position', 0, 0, 100000);

            /**
             * Enable/disable directional light that always points in the direction the user is looking.
             * Removed in X3D V3.3. See NavigationInfo
             * still supported but required changing default to undefined
             * @var {x3dom.fields.SFBool} headlight
             * @memberof x3dom.nodeTypes.GeoViewpoint
             * @initvalue true; undefined since could be already given by NavigationInfo
             * @field x3dom
             * @instance
             */
            this.addField_SFBool(ctx, 'headlight', undefined);

            /**
             * Specifies the navigation type.
             * Removed in X3D V3.3. See NavigationInfo
             * still supported but required changing default to undefined
             * @var {x3dom.fields.MFString} navType
             * @memberof x3dom.nodeTypes.GeoViewpoint
             * @initvalue ['EXAMINE']; undefined since could be already given by NavigationInfo
             * @field x3dom
             * @instance
             */
            this.addField_MFString(ctx, 'navType', undefined);

            /**
             * The speedFactor field of the GeoViewpoint node is used as a multiplier to the elevation-based velocity that the node sets internally; i.e., this is a relative value and not an absolute speed as is the case for the NavigationInfo node.
             * @var {x3dom.fields.SFFloat} speedFactor
             * @range [0, inf]
             * @memberof x3dom.nodeTypes.GeoViewpoint
             * @initvalue 1.0
             * @field x3d
             * @instance
             */
            this.addField_SFFloat(ctx, 'speedFactor', 1.0);
            
            /**
             * Enable/disable elevation scaled speed for automatic adjustment of user movement as recommended in spec. 
             * Custom field to allow disabling for performance
             * @var {x3dom.fields.SFBool} elevationScaling
             * @memberof x3dom.nodeTypes.GeoViewpoint
             * @initvalue true
             * @field x3dom
             * @instance
             */
            this.addField_SFBool(ctx, 'elevationScaling', true);

            /**
             * The geoOrigin field is used to specify a local coordinate frame for extended precision.
             * @var {x3dom.fields.SFNode} geoOrigin
             * @memberof x3dom.nodeTypes.GeoViewpoint
             * @initvalue x3dom.nodeTypes.X3DViewpointNode
             * @field x3d
             * @instance
             */
            this.addField_SFNode('geoOrigin', x3dom.nodeTypes.GeoOrigin);
            
            // save centerOfRotation field for reset
            this._geoCenterOfRotation = this._vf.centerOfRotation ;

        },
        {
            //overwrite activate function to save initial speed
            //from X3DViewpoint.js
            activate: function (prev) {
                var viewarea = this._nameSpace.doc._viewarea;
                
                if (prev) {
                    viewarea.animateTo(this, prev._autoGen ? null : prev);
                }
                viewarea._needNavigationMatrixUpdate = true;
                
                x3dom.nodeTypes.X3DBindableNode.prototype.activate.call(this, prev);
                
                var navi = viewarea._scene.getNavigationInfo();
                this._initSpeed = navi._vf.speed;
                this._examineSpeed = navi._vf.speed;
                this._lastSpeed = navi._vf.speed;
                this._userSpeedFactor = 1.0;
                this._lastNavType = navi.getType();
                  x3dom.debug.logInfo("initial navigation speed: " + this._initSpeed);
                  x3dom.debug.logInfo(this._xmlNode.hasAttribute('headlight'));
                //set headlight and navType here if they are given (but removed from spec.)
                //is there a way to check if fields are given in the document ? (dom has default values if not given)
                if (this._vf.headlight !== undefined) {navi._vf.headlight = this._vf.headlight;}
                if (this._vf.navType !== undefined) {navi._vf.navType = this._vf.navType;}
                
            },
            
            //overwrite deactivate to restore initial speed
            deactivate: function (prev) {
                var viewarea = this._nameSpace.doc._viewarea;
                var navi = viewarea._scene.getNavigationInfo();
                //retain examine mode speed modifications
                navi._vf.speed = this._examineSpeed;
                x3dom.debug.logInfo("restored navigation speed: " + this._examineSpeed);
                x3dom.nodeTypes.X3DBindableNode.prototype.deactivate.call(this, prev);
                //somehow this.getViewMatrix here gets called one more time after deactivate and resets speed, check there
            },
            
            nodeChanged: function() {
                //this is otherwise in X3DBindableNode but function overwritten here
                this._stack = this._nameSpace.doc._bindableBag.addBindable(this);
                
                //for local use
                this._geoOrigin = this._cf.geoOrigin;
                this._geoSystem = this._vf.geoSystem;
                this._position = this._vf.position;
                this._orientation = this._vf.orientation;
                
                //needs to be here because of GeoOrigin subnode
                this._viewMatrix = this.getInitViewMatrix(this._orientation, this._geoSystem, this._geoOrigin, this._position);

                // also transform centerOfRotation for initial view                
                this._vf.centerOfRotation = this.getGeoCenterOfRotation(this._geoSystem, this._geoOrigin, this._geoCenterOfRotation);
                
                // borrowed from Viewpoint.js
            
                this._projMatrix = null;
                this._lastAspect = 1.0;
 
                // z-ratio: a value around 5000 would be better...
                this._zRatio = 10000;
                // set to -1 to trigger automatic setting since fields do not exist
                this._zNear = -1;
                this._zFar = -1;

                // special stuff...
                this._imgPlaneHeightAtDistOne = 2.0 * Math.tan(this._vf.fieldOfView / 2.0);
                
            },
            // all borrowed from Viewpoint.js
            fieldChanged: function (fieldName) {
                if (fieldName == "position" || fieldName == "orientation") {
                    this.resetView();
                }
                else if (fieldName == "fieldOfView" ||
                    fieldName == "zNear" || fieldName == "zFar") {
                    this._projMatrix = null;   // only trigger refresh
                    // set to -1 to trigger automatic setting since fields do not exist
                    this._zNear = -1;
                    this._zFar = -1;
                    this._imgPlaneHeightAtDistOne = 2.0 * Math.tan(this._vf.fieldOfView / 2.0);
                }
                else if (fieldName.indexOf("bind") >= 0) {
                    // FIXME; call parent.fieldChanged();
                    this.bind(this._vf.bind);
                }
            },

            setProjectionMatrix: function(matrix)
            {
                this._projMatrix = matrix;
            },

            getCenterOfRotation: function() {
                // is already transformed to GC
                return this._vf.centerOfRotation;
            },
            
            getGeoCenterOfRotation: function(geoSystem, geoOrigin, geoCenterOfRotation) {
                var coords = new x3dom.fields.MFVec3f();
                coords.push(geoCenterOfRotation);
                var transformed = x3dom.nodeTypes.GeoCoordinate.prototype.GEOtoX3D(geoSystem, geoOrigin, coords);
                return transformed[0];
            },
            
            isExamineMode: function(navType) {
                return (navType == 'examine' || navType == 'turntable' || navType == 'lookaround' || navType == 'lookat');
            },
            
            getViewMatrix: function() {
                //called a lot from viewarea.js; (ab)use for updating elevation scaled speed
                //do only if elevationScaling is enabled
                //skip frames for performance ? do only every 0.1s or so ?
                //gets called once even after being deactivated
                if (this._vf.isActive && this._vf.elevationScaling) {
                    var viewarea = this._nameSpace.doc._viewarea;
                    var navi = viewarea._scene.getNavigationInfo();
                    var navType = navi.getType();
                    //manage examine mode: do not use elevation scaled speed and keep own speed
                    if (this.isExamineMode(navType)) {
                        if (!this.isExamineMode(this._lastNavType)) {
                            x3dom.debug.logError("examine mode speed: " + this._examineSpeed);
                            navi._vf.speed = this._examineSpeed;
                        }
                        this._lastNavType = navType;
                    }
                    else {
                        if (this.isExamineMode(this._lastNavType)) {
                            this._examineSpeed = navi._vf.speed;
                            x3dom.debug.logError("returned from examine: " + this._lastSpeed);
                            navi._vf.speed = this._lastSpeed;
                        }
                        this._lastNavType = navType;
                        //check if speed was modified interactively
                        if (navi._vf.speed != this._lastSpeed) {
                            this._userSpeedFactor *= navi._vf.speed / this._lastSpeed;
                            x3dom.debug.logError("interactive speed factor changed: " + this._userSpeedFactor);
                        }
                        // get elevation above ground
                        // current position in x3d 
                        // borrowed from webgl_gfx.js
                        var viewtrafo = viewarea._scene.getViewpoint().getCurrentTransform();
                        viewtrafo = viewtrafo.inverse().mult(this._viewMatrix);
                        var position = viewtrafo.inverse().e3();
                        
                        ////TODO: deal with geoOrigin here since below only valid for GC
                        ////need inverse geoOrigin; add back offset but how to do rotateYUp: use inverse matrix ?
                        ////eg. first rotate back, then translate back
                        
                        var geoOrigin = this._geoOrigin;
                        var geoSystem = this._geoSystem;
                        var positionGC = position;
                        
                        if (geoOrigin.node) {
                            var origin = x3dom.nodeTypes.GeoCoordinate.prototype.OriginToGC(geoOrigin);
                            if(geoOrigin.node._vf.rotateYUp) {
                                //rotation is GeoLocation rotation
                                var rotmat = x3dom.nodeTypes.GeoLocation.prototype.getGeoRotMat(geoSystem, origin);
                                positionGC = rotmat.multMatrixPnt(position);
                                //return rotMat.inverse().mult(x3dom.fields.SFMatrix4f.translation(transformed.subtract(origin)).mult(rotMat));
                            }
                            //then translate; account for geoOrigin by adding GeoOrigin
                            positionGC = positionGC.add(origin);
                        }
                        
                        // x3dom.debug.logInfo("viewpoint position " + positionGC);
                        
                        var coords = new x3dom.fields.MFVec3f();
                        coords.push(positionGC);
                        //could be a bit optimized since geoSystem does not change
                        var positionGD = x3dom.nodeTypes.GeoCoordinate.prototype.GCtoGD(geoSystem, coords)[0];
                        var elevation = positionGD.z;
                          // x3dom.debug.logInfo("Geoelevation is " + elevation);
                        // at 10m above ground a speed of 1 sounds about right; is positive
                        navi._vf.speed = Math.abs(elevation/10.0) * this._vf.speedFactor * this._userSpeedFactor;
                        this._lastSpeed = navi._vf.speed;
                          // x3dom.debug.logInfo("Changed navigation speed to " + navi._vf.speed + "; ground position at: " + positionGD.y + ", " + positionGD.x);
                    }
                }
                return this._viewMatrix;
            },
            
            getInitViewMatrix: function(orientation, geoSystem, geoOrigin, position) {
                // orientation needs to rotated as in GeoLocation node
                var coords = new x3dom.fields.MFVec3f();
                coords.push(position);
                var positionGC = x3dom.nodeTypes.GeoCoordinate.prototype.GEOtoGC(geoSystem, geoOrigin, coords)[0];
                  //x3dom.debug.logInfo("GEOVIEWPOINT at GC: " + positionGC);
                var orientMatrix = orientation.toMatrix();
                var rotMat = x3dom.nodeTypes.GeoLocation.prototype.getGeoRotMat(geoSystem, positionGC);
                var rotOrient = rotMat.mult(orientMatrix);
                //revert for rotateYUp
                if(geoOrigin.node) {
                    if(geoOrigin.node._vf.rotateYUp) {
                        var origin = x3dom.nodeTypes.GeoCoordinate.prototype.OriginToGC(geoOrigin);
                        var rotMatOrigin = x3dom.nodeTypes.GeoLocation.prototype.getGeoRotMat(geoSystem, origin);
                        rotOrient = rotMatOrigin.inverse().mult(rotOrient);
                    }
                }
                var positionX3D = x3dom.nodeTypes.GeoCoordinate.prototype.GEOtoX3D(geoSystem, geoOrigin, coords)[0];
                  //x3dom.debug.logInfo("GEOVIEWPOINT at X3D: " + positionX3D);
                return x3dom.fields.SFMatrix4f.translation(positionX3D).mult(rotOrient).inverse();
            },

            getFieldOfView: function() {
                return this._vf.fieldOfView;
            },

            resetView: function() {
                this._viewMatrix = this.getInitViewMatrix(this._vf.orientation, this._vf.geoSystem, this._cf.geoOrigin, this._vf.position);
                //also reset center of Rotation ? Not done for regular viewpoint; would need to save original
                this._vf.centerOfRotation = this.getGeoCenterOfRotation(this._vf.geoSystem, this._cf.geoOrigin, this._geoCenterOfRotation);
            },

            getNear: function() {
                return this._zNear;
            },

            getFar: function() {
                return this._zFar;
            },

            getImgPlaneHeightAtDistOne: function() {
                return this._imgPlaneHeightAtDistOne;
            },

            getProjectionMatrix: function(aspect)
            {
                var fovy = this._vf.fieldOfView;
                // set to -1 to trigger automatic setting since fields do not exist
                var zfar = -1;
                var znear = -1;

                if (znear <= 0 || zfar <= 0)
                {
                    var nearScale = 0.8, farScale = 1.2;
                    var viewarea = this._nameSpace.doc._viewarea;
                    var scene = viewarea._scene;

                    // Doesn't work if called e.g. from RenderedTexture with different sub-scene
                    var min = x3dom.fields.SFVec3f.copy(scene._lastMin);
                    var max = x3dom.fields.SFVec3f.copy(scene._lastMax);

                    var dia = max.subtract(min);
                    var sRad = dia.length() / 2;

                    var mat = viewarea.getViewMatrix().inverse();
                    var vp = mat.e3();

                    // account for scales around the viewpoint
                    var translation = new x3dom.fields.SFVec3f(0,0,0),
                        scaleFactor = new x3dom.fields.SFVec3f(1,1,1);
                    var rotation = new x3dom.fields.Quaternion(0,0,1,0),
                        scaleOrientation = new x3dom.fields.Quaternion(0,0,1,0);

                    // unfortunately, decompose is a rather expensive operation
                    mat.getTransform(translation, rotation, scaleFactor, scaleOrientation);

                    var minScal = scaleFactor.x, maxScal = scaleFactor.x;

                    if (maxScal < scaleFactor.y) maxScal = scaleFactor.y;
                    if (minScal > scaleFactor.y) minScal = scaleFactor.y;
                    if (maxScal < scaleFactor.z) maxScal = scaleFactor.z;
                    if (minScal > scaleFactor.z) minScal = scaleFactor.z;

                    if (maxScal > 1)
                        nearScale /= maxScal;
                    else if (minScal > x3dom.fields.Eps && minScal < 1)
                        farScale /= minScal;
                    // near/far scale adaption done

                    var sCenter = min.add(dia.multiply(0.5));
                    var vDist = (vp.subtract(sCenter)).length();

                    if (sRad) {
                        if (vDist > sRad)
                            znear = (vDist - sRad) * nearScale;  // Camera outside scene
                        else
                            znear = 0;                           // Camera inside scene

                        zfar = (vDist + sRad) * farScale;
                    }
                    else {
                        znear = 0.1;
                        zfar = 100000;
                    }

                    var zNearLimit = zfar / this._zRatio;
                    znear = Math.max(znear, Math.max(x3dom.fields.Eps, zNearLimit));
                    //hm, fields do not exist, becomes non-sensical
                    //if (zfar > this._vf.zNear && this._vf.zNear > 0)
                    if (zfar > -1 && -1 > 0)
                        //znear = this._vf.zNear;
                        znear = -1;
                    //if (this._vf.zFar > znear)
                    if (-1 > znear)
                        //zfar = this._vf.zFar;
                        zfar = -1;
                    if (zfar <= znear)
                        zfar = znear + 1;
                    //x3dom.debug.logInfo("near: " + znear + " -> far:" + zfar);
                }

                if (this._projMatrix == null)
                {
                    this._projMatrix = x3dom.fields.SFMatrix4f.perspective(fovy, aspect, znear, zfar);
                }
                else if (this._zNear != znear || this._zFar != zfar)
                {
                    var div = znear - zfar;
                    this._projMatrix._22 = (znear + zfar) / div;
                    this._projMatrix._23 = 2 * znear * zfar / div;
                }
                else if (this._lastAspect != aspect)
                {
                    this._projMatrix._00 = (1 / Math.tan(fovy / 2)) / aspect;
                    this._lastAspect = aspect;
                }

                // also needed for being able to ask for near and far
                this._zNear = znear;
                this._zFar = zfar;

                return this._projMatrix;
            }
        }
    )
);