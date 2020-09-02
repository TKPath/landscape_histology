// to create classified image using trained classifier and determine class % within tissue area defined by simple tissue threshold
// change all names of classifier

//Tim Kendall, University of Edinburgh
//Tim.Kendall@ed.ac.uk

// August 2020

import static qupath.lib.gui.scripting.QPEx.*

setImageType('BRIGHTFIELD_OTHER');
setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049 ", "Stain 2" : "DAB", "Values 2" : "0.26917 0.56824 0.77759 ", "Background" : " 255 255 255 "}');
resetSelection();

// apply tissue threshold to create an annotation, fill holes
createAnnotationsFromPixelClassifier("tissue_whole", 1000000.0, 10000000.0, "SELECT_NEW")

// export the labelled annotation that measurements are restricuerd to
def imageData = getCurrentImageData()

// Define output path (relative to project)
def outputDir = buildFilePath(PROJECT_BASE_DIR, 'parent_annotation')
mkdirs(outputDir)
def name_full = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())
def path_full = buildFilePath(outputDir, name_full + "_parent.png")

// Define how much to downsample during export (may be required for large images)
double downsample = 8

// Create an ImageServer where the pixels are derived from annotations
def labelServer = new LabeledImageServer.Builder(imageData)
  .backgroundLabel(0, ColorTools.WHITE) // Specify background label (usually 0 or 255)
  .downsample(downsample)    // Choose server resolution; this should match the resolution at which tiles are exported
  .addLabel('tissue', 1)      // Choose output labels (the order matters!)
  .multichannelOutput(false) // If true, each label refers to the channel of a multichannel binary image (required for multiclass probability)
  .build()

// Write the image
writeImage(labelServer, path_full)
print 'Parent annotation image exported to ' + path_full

// measure based on thyroid classifier only within the tissue annotation
addPixelClassifierMeasurements("thyroid_whole2", "thyroid_whole2")

// Save the classified image
def server = getCurrentServer()
def name = server.getMetadata().getName()
def name_short = name.reverse().drop(4).reverse() + "_classified.tif"
def classified_path = buildFilePath(PROJECT_BASE_DIR, 'classified_images')
mkdirs(classified_path)
pathClassified = buildFilePath(classified_path, name_short)
writePredictionImage('thyroid_whole2', pathClassified)
print 'Image exported to ' + pathClassified

// Export the results - Counts
def name_count = name.reverse().drop(4).reverse() + ".csv"
def count_path = buildFilePath(PROJECT_BASE_DIR, 'counts')
mkdirs(count_path)
pathCount = buildFilePath(count_path, name_count)
saveAnnotationMeasurements(pathCount)
print 'Results exported to ' + pathCount
