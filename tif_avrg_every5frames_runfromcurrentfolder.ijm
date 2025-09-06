// Take tif files, make avrg z projection every 5 frames,concatenate into a time series

//// 
rootdir = getDirectory("Choose directory containing input files")
input = rootdir+ "TIF/"; 
output = rootdir+"TIF_avrg5frame/"

File.makeDirectory(output);
function action(input, output, filename) {
		// change the total number of frames for the down-sampled images
		//nframe_out = 15;// 75 frames at average 5 per frame (10 minute video, 8s per frame)
		//nframe_out = 24;// 120 frames at average 5 per frame (16 minute video, 8s per frame)
		nframe_out = 36;// 180 frames at average 5 per frame (15 minute video, 5s per frame)
		open(input + filename);
		selectWindow(filename);
		run("Duplicate...", "duplicate range=1-5 use");
		run("Z Project...", "projection=[Average Intensity]");
		selectWindow(filename);
		run("Duplicate...", "duplicate range=6-10");
		run("Z Project...", "projection=[Average Intensity]");
		title_out = substring(filename,0,lengthOf(filename) -4) +"_AVG";
		title1 = "AVG_"+ substring(filename,0,lengthOf(filename)-4) +"-1.tif";
		title2 = "AVG_"+ substring(filename,0,lengthOf(filename) -4) +"-2.tif";
		run("Concatenate...", "  title=" + title_out +" open image1="+title1+" image2="+title2+" image3=[-- None --]");
		selectWindow(substring(title1,4,lengthOf(title1)));
		close();
		selectWindow(substring(title2,4,lengthOf(title2)));
		close();

		for (j = 2; j < nframe_out; j++){
        selectWindow(filename);
        start_frame = 1+j*5;
        end_frame = (j+1)*5;
        run("Duplicate...", "duplicate range="+start_frame + "-"+ end_frame); 
        run("Z Project...", "projection=[Average Intensity]");
        
        run("Concatenate...", "  title="+title_out+" open image1="+title_out+" image2="+title1+" image3=[-- None --]");
		selectWindow(substring(title1,4,lengthOf(title1)));
		close();
		}	
		selectWindow(title_out);
		saveAs("Tiff", output + title_out+".tif");
		close();
		close();
}

list = getFileList(input);

for (i = 0; i < list.length; i++){
		if (endsWith(list[i], ".tif")){
        	action(input, output, list[i]);
		}
}

