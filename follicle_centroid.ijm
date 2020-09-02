//Auto follicle centroid determination
//Saves outlines and csv of xy centroids
//Tim Kendall, University of Edinburgh
//Tim.Kendall@ed.ac.uk

// August 2020

dir = getDirectory("Choose a Directory ");
count = 1; 
run("Set Measurements...", "centroid redirect=None decimal=3");
follicle_centroid(dir)


function follicle_centroid(dir) {
	list = getFileList(dir);      	
	for (i=0; i<list.length; i++) {
      		open(dir+list[i]);
      		b=getTitle();
         	print(b);
         	setThreshold(2, 2); // select the colloid mask
         	setOption("BlackBackground", false);
			run("Convert to Mask");
			selectWindow(b);
			run("Fill Holes");
			run("Analyze Particles...", "size=100-Infinity circularity=0.10-1.00 show=Outlines display");
			saveAs("Results", dir+File.separator+b+"_centroids.csv");
			selectWindow("Results");
			run("Close");
			selectWindow("Drawing of "+b);
			saveAs("PNG", dir+File.separator+b+"_outlines");
			close();
			selectWindow(b);
			close();
      	}
}