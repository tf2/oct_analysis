get_max_peak_trove <- function(up_down, l) {
	ds = vector()
	for(x in 2:(length(up_down)-1)) {
		if(l$y[up_down[x-1]]> l$y[up_down[x]] & l$y[up_down[x+1]]> l$y[up_down[x]])
			ds[x] = abs(l$y[up_down[x]] - l$y[up_down[x-1]]) + abs(l$y[up_down[x]] - l$y[up_down[x+1]])
	}
	ds[1] = 0; ds = c(ds, 0)
	params = NULL
	params$midpoint = up_down[which.max(ds)]
	params$l_shoulder = up_down[which.max(ds)-1]
	params$r_shoulder = up_down[which.max(ds)+1]
return(params)
}

check_up_and_down <- function(l) {
	check = 1
	up_down = NULL
	for(x in 2:length(l$x)) {
		if(check==1 & l$y[x]<l$y[x-1]) {
			check = 2
			up_down = c(up_down, l$x[x])
		}
		if(check==2 & l$y[x]>l$y[x-1]) {
			check = 1
			up_down = c(up_down, l$x[x])
		}
	}
return(up_down)
}

get_fovea_position <- function(img, fac=0.2) {
	r = NULL
	for(x in 1:dim(img)[3]) {
		a = apply(img[,,x], 2, mean)
		r = rbind(r, a)
	}
	l = lowess(apply(r, 2, mean), f=fac)
	up_down = check_up_and_down(l)
	fpos = get_max_peak_trove(up_down, l)$midpoint
return(fpos)
}

get_fovea_slices<- function(img, fpos, fac=0.2) {
	r = NULL
	for(x in 1:dim(img)[1]) {
		a = apply(img[x,(fpos-5):(fpos+5),], 2, mean)
		r = rbind(r, a)
	}
	l = lowess(apply(r, 2, mean), f=fac)
	up_down = check_up_and_down(l)
	fpos = get_max_peak_trove(up_down, l)
return(fpos)
}

fit_x_fovea <- function(fovea_slice, fac=0.2) {
	x_dist = apply(fovea_slice, 1, mean)
	x_dist = x_dist-median(x_dist)
	l = lowess(x_dist, f=fac)
	up_down = check_up_and_down(l)
	fovea_params = get_max_peak_trove(up_down, l)
return(fovea_params)
}

fit_y_fovea <- function(fovea_slice, x_fovea_params, fac=0.1) {
	x_dist = apply(fovea_slice[x_fovea_params$l_shoulder:x_fovea_params$r_shoulder,], 2, mean)
	x_dist = x_dist-median(x_dist)
	x_dist[ which.min(x_dist):1] = NA
	params = NULL
	params$dip_pos = which.max(x_dist)
	x_dist[ params$dip_pos:1] = NA
	l = lowess(x_dist, f=fac)
	params$shoulder_r_pos = min(which(l$y<0))
	left = fovea_slice[x_fovea_params$l_shoulder:x_fovea_params$midpoint,params$dip_pos:dim(fovea_slice)[2]]
	lef = apply(left, 2, mean)
	lef = lef-median(lef)
	mp = check_up_and_down(lowess(cumsum(lef), f=0.01))[1]

	x1_adjust = check_up_and_down(lowess(cumsum(left[,mp]), f=0.01))
	params$left_position = x_fovea_params$l_shoulder
	if(!is.null(x1_adjust)) {
		params$left_position = params$left_position + max(x1_adjust)
	}
	x_fovea_params$l_shoulder = params$left_position
	#cut = quantile(abs(lef), probs=0.68)
	#mp = max(which(lef  > median(lef)+(cut*3)))
	#mp =  max(which(lowess (lef)$y  > median(lef)+(cut*3)))
	#l = lowess(lef, f=0.5)
	#mp = which.max(lm(l$y~l$x)$residuals)
	params$left_top_sholder = params$dip_pos+mp
	right = fovea_slice[x_fovea_params$midpoint:x_fovea_params$r_shoulder,params$dip_pos:dim(fovea_slice)[2]]
	rig = apply(right, 2, mean)
	rig = rig-median(rig)
	mp = check_up_and_down(lowess(cumsum(rig), f=0.01))[1]

	x2_adjust = check_up_and_down(lowess(cumsum(right[,mp]), f=0.01))
	params$right_position = x_fovea_params$r_shoulder
	if(!is.null(x2_adjust)) {
		params$right_position = params$right_position - min(x2_adjust)
	}
	x_fovea_params$r_shoulder = params$right_position
	#mp=  max(which(lowess (right)$y  > median(rig)+sd(rig)))
	#l = lowess(rig, f=0.5)
	#which.max(lm(l$y ~ poly(l$x,3))$residuals)
	#mp = which.max(lm(l$y~l$x)$residuals)
	params$right_top_sholder = params$dip_pos+mp
return(params)
}

fit_sholders <- function(fovea_slice, x_fovea_params, y_fovea_params) {
	#max_y = max(c(y_fovea_params$left_top_sholder, y_fovea_params$right_top_sholder))
	#zoom_fovea = fovea_slice[x_fovea_params$l_shoulder:x_fovea_params$r_shoulder,y_fovea_params$dip_pos:max_y]
	zoom_fovea = fovea_slice[x_fovea_params$l_shoulder:x_fovea_params$r_shoulder,y_fovea_params$dip_pos:y_fovea_params$shoulder_r_pos]
	y_dist = apply(zoom_fovea, 1, mean)
	l = lowess(y_dist, f=0.1)
	#image(zoom_fovea)
	#matplot(l$x/dim(zoom_fovea)[1], l$y/dim(zoom_fovea)[2], type='l', col="blue", add=T)

	#image(fovea_slice)
	scales = dim(zoom_fovea)/dim(fovea_slice)
	x_add = x_fovea_params$l_shoulder/dim(fovea_slice)[1]
	y_add = y_fovea_params$dip_pos/ dim(fovea_slice)[2]
	#matplot(((l$x/dim(zoom_fovea)[1])*scales[1])+x_add, ((l$y/dim(zoom_fovea)[2])*scales[2])+y_add, type='l', col="blue", add=T)
return(list("x" = ((l$x/dim(zoom_fovea)[1])*scales[1])+x_add, "y" = ((l$y/dim(zoom_fovea)[2])*scales[2])+y_add))
}

plot_fitted_slice <- function(fovea_slice, x_fovea_params, y_fovea_params, dip) {
	p_ind = dim(fovea_slice)
	m1 = x_fovea_params$midpoint/p_ind[1]
	x1 = x_fovea_params$l_shoulder/p_ind[1]
	x2 = x_fovea_params$r_shoulder/p_ind[1]
	
	# this for dynamic x points
	#x1 = y_fovea_params$left_position/p_ind[1]
	#x2 = y_fovea_params$right_position/p_ind[1]
	
	y1 = y_fovea_params$dip_pos/p_ind[2]
	y2 = y_fovea_params$shoulder_r_pos/p_ind[2]
	ly1 = y_fovea_params$left_top_sholder/p_ind[2]
	ry1 = y_fovea_params$right_top_sholder/p_ind[2]
	ly2 = y1
	ry2 = y1
	coor = c(ly1, ry1)
	if(which.min(coor)==1) {
		ly2 = ly2 - (ry1-ly1)
		ry2 = ry2 + (ry1-ly1)
	} else {
		ly2 = ly2 + (ry1-ly1)
		ry2 = ry2 - (ry1-ly1)
	}

	image(fovea_slice, breaks=c(0,50,200), col=c("black", "green"))
	#segments(x1,ly2,x2,ry2, col="blue")
	segments(x1,y1,x2,y1, col="blue")
	segments(x1,ly1,x2,ry1, col="blue")
	segments(x1,y1,x1,ly1, col="blue")
	segments(x2,y1,x2,ry1, col="blue")
	segments(m1,y1,m1,mean(c(ly1, ry1)), col="red")

	co = NULL
	if(length(dip$y[dip$y<mean(c(ly1, ry1))])/length(dip$y)>0.2) {
		co = "red"
		matplot(dip$x, dip$y, lwd=2, type='l', col=co, add=T)
	}
	if(length(dip$y[dip$y<mean(c(ly1, ry1))])/length(dip$y)>0.5) {
		co = "blue"
		matplot(dip$x, dip$y, lwd=2, type='l', col=co, add=T)
	}
}

library(oro.dicom)
source("dicom.R")
smooth_factor = 0.2
image_file = "../data/656057.dcm"
image_data = dicomInfo(image_file)

fpos = get_fovea_position(image_data$img)
slices = get_fovea_slices(image_data$img, fpos[1])
slice = slices$midpoint

fovea_slice = t(image_data$img[,,slice])
x_fovea_params = fit_x_fovea(fovea_slice)
y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)
dip = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)

plot_fitted_slice(fovea_slice, x_fovea_params, y_fovea_params, dip)





###############################
### Fit across the foveal dip - single dataset 
cut = 0.5
pdf("try_fit_center_test.pdf", onefile=T)
for(slice in slices$l_shoulder:slices$r_shoulder) {
	fovea_slice = t(image_data$img[,,slice])
	y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)
	dip = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)
	plot_fitted_slice(fovea_slice, x_fovea_params, y_fovea_params, dip)
}
dev.off()


### All max foveal dips - across files
smooth_factor = 0.2
pdf("all_foveas.pdf", onefile=T)
direct = "../data/"
files = dir(direct)
for(x in 1:length(files)) {
	image_file = paste(direct, "/", files[x], sep="")
	image_data = dicomInfo(image_file)
	
	fpos = get_fovea_position(image_data$img)
	slices = get_fovea_slices(image_data$img, fpos[1])
	slice = slices$midpoint

	fovea_slice = t(image_data$img[,,slice])
	x_fovea_params = fit_x_fovea(fovea_slice)
	y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)
	dip = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)

	plot_fitted_slice(fovea_slice, x_fovea_params, y_fovea_params, dip)
}
dev.off()

#####


decend_forward <- function(co, span=3) {
	for(x in 1:(length(co)-3) ) {
		pin = 0
		for(y in 0:span) {
			if(co[x+y]=="red") {
				pin = pin+1
			}
		}
		if(pin!=span) {
			co[x] = "black"
		}
	}
return(co)
}

decend_backward <- function(co) {
	for(x in (length(co):3) ) {
		if (co[x]=="red" & co[x-1]=="red" & co[x-2]=="red") {
			co[x] = "red"
			co[x-1] = "red"
			co[x-2] = "red"
		} else {
			co[x] = "black"
		}
	}
return(co)
}



registar <- function(img) {
	r = NULL
	par(mfrow=c(1,2))
	fac = ceiling(dim(img)[3]/4)
	for(x in 1:dim(img)[3]) {
		a = apply(img[,,x], 1, median)
		bands = img[which(a> median(a)+sd(a)),,x]
		f = a[order(a, decreasing=T)]
		#bands1 = img[a>=f[fac],,x]
		#bands1 = bands1[1:fac,]
		image(t(bands))
		image(t(img[,,x]))
		scan("")
	}
}

denoise <- function(mat, span=2) {
	for(x in 1:length(mat[,1])) {
		m = mat[x,]
		s = seq(1, length(m), by=span)
		for(y in 1:(length(s)-1)) {
			if ( sum(!is.na(m[s[y]:(s[y]+span)])) != span ) {
				m[s[y]:(s[y]+span)] = NA
			}
		}
		mat[x,] = m
	}
}

fit_bands <- function(img) {
	qunats = vector()
	img = t(image_data$img[,,slice])
	ma = 0
	mat = img
	for(x in 1:dim(img)[1]) {
		l  = lowess(img[x,], f=0.05)
		if(max(img[x,])>ma) {
			ma = max(img[x,])
		}
		l$y = l$y-median(l$y)
		#qunats[x] = quantile(abs(l$y), probs=0.68)*5
		cut = sd(l$y)*3
		m = mat[x,]
		m[m<cut]  = NA
		m[m>cut]  = 100
		mat[x,] = m
		qunats[x] = cut
	}
	

	image(img, breaks=c(0,mean(qunats), ma), col=c("black", "green"))
	image(mat, add=T, col="red")
}
