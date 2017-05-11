#' Transforming your WGS84 data to TWD97 data
#'
#' This function allows you to transform your WGS84 location data to TWD97 data and create a new data frame.
#' WARNING: ONLY FOR LOCATION OF TAIWAN!!!!!!
#' @param inputLat A numeric vector, WGS84 latitude of yours.
#' @param inputLng A numeric vector, WGS84 longitude of yours.
#' @keywords WGS84toTWD97
#' @export
#' @examples
#' newTWD97data <- WGS84toTWD97(inputLat = 24.321, inputLng = 121.234)

WGS84toTWD97 <- function(inputLat = NULL, inputLng = NULL){
  # inputLat
  # inputLng 
  # ============= Calculate Section =============
  #  TWD97：
  # 台灣使用的是2度分帶，以東經121度為中央經線。
  
  # base constant for our function
  
  a <- 6378137.0 # Equatorial Radius (地球赤道半徑) = 6378137.0 M
  b <- 6356752.314245 # Polar Radius (地球兩極半徑) = 6356752.314245 M
  if(inputLng<120){
    lng0 <- 119   # Penghu, Kinmen and Matsu central longitude (澎湖 金門 馬祖中心經度) = 119 degree
  }else{
    lng0 <- 121 # Taiwan central longitude (台灣中心經度) = 121 degree
  }
  
  k0 <- 0.9999 # scaling size (縮放比例)
  
  dx <- 250000 # 橫座標平移量
  dy <- 0 # 縱坐標平移量
  
  # 參數：
  e <- sqrt(1 - ((b*b)/(a*a)))  # 橢圓面積離心率 1<e<0
  
  f <- as.numeric(a-b)/a # 扁平率
  n <- (a-b)/(a+b)
  # n =parseFloat( (parseFloat(a)-parseFloat(b))/(parseFloat(a)+parseFloat(b)))
  
  #   ~~~~~~~Kruger Series~~~~~~~
  
  # AA = =(a/( 1 + n ) )*(1 + (1/4)* n**2 + (1/64)*n**4 + (1/256)*n**6 + (25/16384)*n**8 + (49/65536)*n**10)
  AA <- (a/( 1 + n ) )*(1 + (1/4)* n**2 + (1/64)*n**4 + (1/256)*n**6 + (25/16384)*n**8 + (49/65536)*n**10)
  
  alpha1 <- (1/2)*n - (2/3)*n**2 + (5/16)*n**3 + (41/180)*n**4 - (127/288)*n**5 + (7891/37800)*n**6 + (72161/387072)*n**7 - (18975107/50803200)*n**8 + (60193001/290304000)*n**9 + (134592031/1026432000)*n**10
  alpha2 <- (13/48)*n**2 - (3/5)*n**3 + (557/1440)*n**4 + (281/630)*n**5 - (1983433/1935360)*n**6 + (13769/28800)*n**7 + (148003883/174182400)*n**8 - (705286231/465696000)*n**9 + (1703267974087/3218890752000)*n**10
  alpha3 <- (61/240)*n**3 - (103/140)*n**4 + (15061/26880)*n**5 + (167603/181440)*n**6 - (67102379/29030400)*n**7 + (79682431/79833600)*n**8 + (6304945039/2128896000)*n**9 -  (6601904925257/1307674368000)*n**10
  alpha4 <- (49561/161280)*n**4 - (179/168)*n**5 + (6601661/7257600)*n**6 + (97445/49896)*n**7 - (40176129013/7664025600)*n**8 + (138471097/66528000)*n**9 + (48087451385201/5230697472000)*n**10
  alpha5 <- (34729/80640)*n**5 - (3418889/1995840)*n**6 + (14644087/9123840)*n**7 +   (2605413599/622702080)*n**8 - (31015475399/2583060480)*n**9 +  (5820486440369/1307674368000)*n**10
  alpha6 <- (212378941/319334400)*n**6 - (30705481/10378368)*n**7 + (175214326799/58118860800)*n**8 + (870492877/96096000)*n**9 - (1328004581729000/47823519744000)*n**10
  alpha7 <- (1522256789/1383782400)*n**7 - (16759934899/3113510400)*n**8 + (1315149374443/221405184000)*n**9 + (71809987837451/3629463552000)*n**10
  alpha8 <- (1424729850961/743921418240)*n**8 -   (256783708069/25204608000)*n**9 + (2468749292989890/203249958912000)*n**10
  alpha9 <- (21091646195357/6080126976000)*n**9 - (67196182138355800/3379030566912000)*n**10
  alpha10 <- (77911515623232800/12014330904576000)*n**10
  
  # Load Latitude and Longitude
  lat <- inputLat # 浮點數計算
  lng <- inputLng # 浮點數計算
  latr <- lat*pi/180 # 弧度
  Dlng <- lng-lng0  # 經度與中央經線相差值
  Dlngr<- Dlng*pi/180 # 弧度
  # conformal latitude
  confLat <- atan(sinh(asinh(tan(latr))-e*atanh(e*sin(latr))))
  
  # sigma
  sigma = sinh(e*atanh(e*tan(latr)/sqrt(1+tan(latr)**2)))
  
  # tau = tan(lat) , taup = tau' = tan(conLat) 
  tau <- tan(latr)
  taup = tan(confLat)
  # xi = North direction, conformal Xi', 
  # eta =  East direction, conformal eta.
  xip <- atan(taup/cos(Dlngr))
  etap <- asinh(sin(Dlngr)/sqrt(taup*taup+(cos(Dlngr)**2)))
  
  xi <- xip+alpha1*sin(2*xip)*cosh(2*etap)+alpha2*sin(4*xip)*cosh(4*etap)+alpha3*sin(6*xip)*cosh(6*etap)+alpha4*sin(8*xip)*cosh(8*etap)+alpha5*sin(10*xip)*cosh(10*etap)+alpha6*sin(12*xip)*cosh(12*etap)+alpha7*sin(14*xip)*cosh(14*etap)
  
  eta <- etap+alpha1*cos(2*xip)*sinh(2*etap)+alpha2*cos(4*xip)*sinh(4*etap)+alpha3*cos(6*xip)*sinh(6*etap)+alpha4*cos(8*xip)*sinh(8*etap)+alpha5*cos(10*xip)*sinh(10*etap)+alpha6*cos(12*xip)*sinh(12*etap)+alpha6*cos(14*xip)*sinh(14*etap)
  
  easting <- k0*AA*eta
  northing <- k0*AA*xi
  
  realEasting <- dx+easting
  WGS84toTWD97.dataframe <- data.frame(Latitude_TWD97_Y = northing, Longitude_TWD97_X = realEasting)
  
  return(WGS84toTWD97.dataframe)
}

#' Transforming your TWD97 data to WGS84 data
#'
#' This function allows you to transform your TWD97 location data to WGS84 data and create a new data frame.
#' WARNING: ONLY FOR LOCATION OF TAIWAN!!!!!!
#' @param inputTWD97Y A numeric vector, TWD97 Y of yours.
#' @param inputTWD97X A numeric vector, TWD97 X of yours.
#' @param Lng0 A numeric vector. If your want to use this function for Penghu, Kinmen and Matsu, you must input new central longitude for 119.
#' @keywords TWD97toWGS84
#' @export
#' @examples
#' In Taiwan:
#' transTWD97data <- TWD97toWGS84(inputTWD97Y = 2767181, inputTWD97X = 296988.9)
#' 
#' In Penghu, Kinmen and Matsu:
#' transTWD97data <- TWD97toWGS84(inputTWD97Y = 2703564, inputTWD97X = 180914.8, Lng0 = 119)

TWD97toWGS84 <- function(inputTWD97Y = NULL, inputTWD97X = NULL, Lng0 = 121){
  
  # ============= Calculate Section =============
  a <- 6378137.0 #Equatorial Radius (地球赤道半徑) = 6378137.0 M
  b <- 6356752.314245 #Polar Radius (地球兩極半徑) = 6356752.314245 M
  lng0 <- Lng0 #Taiwan central longitude (台灣中心經度) = 121 degree
  k0 <- 0.9999 # scaling size (縮放比例)
  
  dx <- 250000 #橫座標平移量
  dy <- 0 #縱坐標平移量
  
  # 參數：
  e <- sqrt(1 - ((b*b)/(a*a)))  # 橢圓面積離心率 1<e<0
  
  
  # f = (parseFloat(a)-parseFloat(b))/parseFloat(a) #扁平率
  
  n <- (a-b)/(a+b) 
  AA <- (a/( 1 + n ) )*(1 + (1/4)* n**2 + (1/64)*n**4 + (1/256)*n**6 + (25/16384)*n**8 + (49/65536)*n**10)
  
  #  ~~~~~~~Kruger Series~~~~~~~
  beta1 <- (1/2)*n-(2/3)*n**2 + (37/96)*n**3 - (1/360)*n**4 - (81/512)*n**5 + (96199/604800)*n**6 - (5406467/38707200)*n**7 + (7944359/67737600)*n**8 - (7378753979/97542144000)*n**9 + (25123531261/804722688000)*n**10
  beta2 <- (1/48)*n**2 + (1/15)*n**3 - (437/1440)*n**4 + (46/105)*n**5 - (1118711/3870720)*n**6 + (51841/1209600)*n**7 + (24749483/348364800)*n**8 - (115295683/1397088000)*n**9 + (5487737251099/51502252032000)*n**10
  beta3 <- (17/480)*n**3 - (37/840)*n**4 - (209/4480)*n**5 + (5569/90720)*n**6 + (9261899/58060800)*n**7 - (6457463/17740800)*n**8 + (2473691167/9289728000)*n**9 - (852549456029/20922789888000)*n**10
  beta4 <- (4397/161280)*n**4 - (11/504)*n**5 - (830251/7257600)*n**6 + (466511/2494800)*n**7 + (324154477/7664025600)*n**8 - (937932223/3891888000)*n**9 - (89112264211/5230697472000)*n**10
  beta5 <- (4583/161280)*n**5 - (108847/3991680)*n**6 - (8005831/63866880)*n**7 + (22894433/124540416)*n**8 + (112731569449/557941063680)*n**9 - (5391039814733/10461394944000)*n**10
  beta6 <- (20648693/638668800)*n**6 -  (16363163/518918400)*n**7 - (2204645983/12915302400)*n**8 + (4543317553/18162144000)*n**9 + (54894890298749/167382319104000)*n**10
  beta7 <- (219941297/5535129600)*n**7 - (497323811/12454041600)*n**8 - (79431132943/332107776000)*n**9 + (4346429528407/12703122432000)*n**10
  beta8 <- (191773887257/3719607091200)*n**8 -  (17822319343/336825216000)*n**9 - (497155444501631/1422749712384000)*n**10
  beta9 <- (11025641854267/158083301376000)*n**9  - (492293158444691/6758061133824000)*n**10
  beta10 <- (7028504530429620/72085985427456000)*n**10
  
  # xi and eta
  xi <- inputTWD97Y/(k0*AA)
  eta <- (inputTWD97X-dx)/(k0*AA)
  xip <- xi-(beta1*sin(2*xi)*cosh(2*eta)+beta2*sin(4*xi)*cosh(4*eta)+beta3*sin(6*xi)*cosh(6*eta)+beta4*sin(8*xi)*cosh(8*eta)+beta5*sin(10*xi)*cosh(10*eta)+beta6*sin(12*xi)*cosh(12*eta)+beta7*sin(14*xi)*cosh(14*eta))
  etap <- eta-(beta1*cos(2*xi)*sinh(2*eta)+beta2*cos(4*xi)*sinh(4*eta)+beta3*cos(6*xi)*sinh(6*eta)+beta4*cos(8*xi)*sinh(8*eta)+beta5*cos(10*xi)*sinh(10*eta)+beta6*cos(12*xi)*sinh(12*eta)+beta7*cos(14*xi)*sinh(14*eta))
  taup <- sin(xip)/sqrt(sinh(etap)**2+cos(xip)**2)
  
  # Calculate for Lngitude
  lngr <- atan(sinh(etap)/cos(xip))
  lngd <- lngr*180/pi
  resultLng <- lngd+lng0
  
  # Calculate foe Latitude
  sigma0 <- sinh(e*atanh(e*taup/sqrt(1+taup**2)))
  
  f <- taup*sqrt(1+sigma0**2)-sigma0*sqrt(1+taup**2)-taup
  
  dfTauDtau <- (sqrt((1+sigma0**2)*(1+taup**2))-sigma0*taup)*(1-e**2)*sqrt(1+taup**2)/(1+(1-e**2)*taup**2)
  taup1 <- taup-f/dfTauDtau
  sigma1 <- sinh(e*atan(e*taup1/sqrt(1+taup1**2)))
  taupc<- (taup1-f)/(sqrt((1+taup1**2)*(1+taup1**2))-sigma1*taup1)*(1-e**2)*sqrt(1+taup1**2)/(1+(1-e**2)*taup1**2)
  resultLat <- atan(taup1)*180/pi
  
  TWD97toWGS84.dataframe <- data.frame(Latitude = resultLat, Longitude = resultLng)
  return(TWD97toWGS84.dataframe)
}

#' Use TWD97 location data to calculate the distance between two points
#'
#' This function allows you to use TWD97 location data to calculate the distance between two points, create a new vector.
#' WARNING: ONLY FOR LOCATION OF TAIWAN!!!!!!
#' @param origin_TWD97_X A numeric vector, one of the points' TWD97 X data.
#' @param origin_TWD97_Y A numeric vector, one of the points' TWD97 Y data.
#' @param observation_TWD97_X A numeric vector, another of the points' TWD97 X data.
#' @param observation_TWD97_X A numeric vector, another of the points' TWD97 Y data.
#' @keywords distance_TWD97
#' @export
#' @examples
#' distanceNEW <- distance_TWD97(origin_TWD97_X = 296988.9, origin_TWD97_Y = 2767181, observation_TWD97_X = 196302.7, observation_TWD97_Y = 2828646)

distance_TWD97 <- function(origin_TWD97_X = NULL, 
                           origin_TWD97_Y = NULL, 
                           observation_TWD97_X = NULL, 
                           observation_TWD97_Y = NULL) {
  # origin point location data for TWD97 (x,y)
  # observation point location data for TWD97 (x,y)
  distance <- sqrt((origin_TWD97_X-observation_TWD97_X)*(origin_TWD97_X-observation_TWD97_X)+(origin_TWD97_Y-observation_TWD97_Y)*(origin_TWD97_Y-observation_TWD97_Y))/1000
  return(distance)
}

#' Use WGS84 location data to calculate the distance between two points
#'
#' This function allows you to use WGS84 location data to calculate the distance between two points, create a new vector.

#' @param origin_lat A numeric vector, one of the points' latitude data.
#' @param origin_lng A numeric vector, one of the points' longitude data.
#' @param obervation_lat A numeric vector, another of the points' latitude data.
#' @param obervation_lng A numeric vector, another of the points' longitude data.
#' @keywords distance_WGS84
#' @export
#' @examples
#' distanceNEW <- distance_TWD97(origin_lat = 25.01194, origin_lng = 121.4656, obervation_lat = 25.56660, obervation_lng = 120.4656)


distance_WGS84 <- function(origin_lat = NULL, 
                           origin_lng = NULL, 
                           obervation_lat = NULL, 
                           obervation_lng = NULL) {
  # input origin point location 
  origin_point_pi <- c((origin_lat/180)*pi, (origin_lng/180)*pi)
  # input observation point location data
  obervation_lat_pi <- (obervation_lat/180)*pi
  obervation_lng_pi <- (obervation_lng/180)*pi
  D <- 6378100*acos(sin(obervation_lat_pi)*sin(origin_point_pi[[1]])+cos(obervation_lat_pi)*cos(origin_point_pi[[1]])*cos(obervation_lng_pi-origin_point_pi[[2]]))
  distance <- D/1000
  return(distance)
}



