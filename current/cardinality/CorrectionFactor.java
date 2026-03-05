package cardinality;

import java.io.File;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import structures.FloatList;

/**
 * Correction factors for estimators.
 * @author Brian Bushnell, Chloe
 * @date March 5, 2026
 */
public class CorrectionFactor{
	
	/*--------------------------------------------------------------*/
	/*----------------            Main              ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) {
		String fname=args[0];
		int buckets=(args.length>1 ? Parse.parseIntKMG(args[0]) : 2048);
		initialize(fname, buckets);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static synchronized void initialize(String fname, int buckets) {
		correctionFile=fname;
		correctionBuckets=buckets=CardinalityTracker.powerOf2AtLeast(buckets);
		float[][] matrix=loadFile(fname, buckets);
		CF_MATRIX=matrix;
	}
	
	public static synchronized float[][] loadFile(String fname, int buckets){
		String path=(fname.indexOf('?')>=0 ? Data.findPath(fname) : fname);
		File f=new File(path);
		if(!f.canRead() || !f.isFile()) {
			System.err.println("Can't read "+fname);
			return null;
		}
		FloatList[] lists=null;
		FileFormat ff=FileFormat.testInput(fname, null, false);
		ByteFile bf=ByteFile.makeByteFile(ff, buckets);
		LineParser1 lp=new LineParser1('\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			lp.set(line);
			if(lp.startsWith('#')) {//header
			}else {
				if(lists==null) {lists=new FloatList[lp.terms()];}
				for(int i=0; i<lp.terms(); i++) {
					lists[i].add(lp.parseFloat(i));
				}
			}
		}
		return interpolate(lists, buckets);
	}
	
	public static float[][] interpolate(FloatList[] lists, int buckets) {
		float[][] matrix=new float[lists.length][buckets+1];
		if(lists[0].size()==matrix[0].length) {//TODO: easy mode
			
		}else if(lists[0].size()>matrix[0].length) {//TODO: Skip some
			
		}else {//TODO: interpolate using average of the quadratic fits of the pair of 3 closest points
			
		}
		return matrix;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Accessors           ----------------*/
	/*--------------------------------------------------------------*/
	
	public static float getCF(int occupancy, int buckets, int type) {
		float[] array=CF_MATRIX[type];
		if(buckets<=correctionBuckets) {return array[occupancy*(correctionBuckets/buckets)];}
		//TODO: Return the correct CF or interpolate it if there are more buckets than correctionBuckets
		return 1;
	}
	
	public static float getFullCF(int occupancy, int buckets, int type) {
		float[] array=CF_MATRIX[type];
		return array[array.length-1];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static String correctionFile="?cardinalityCorrection.tsv.gz";
	public static int correctionBuckets=2048;
	public static final int OCCUPIED=0, MEAN=1, HMEAN=2, GMEAN=3, HLL=4, LINEAR=5;
	
	
	public static float[][] CF_MATRIX=loadFile(correctionFile, correctionBuckets);
	
}
