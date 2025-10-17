package illumina;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.KillSwitch;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Converts Illumina CBCL files to tab-delimited text format.
 * Reads s.locs for positions, .filter files for pass-filter flags,
 * and .cbcl files for base calls and quality scores.
 *
 * Usage: Cbcl2Text runfolder=<path> out=<file> lane=<int> [tiles=<list>]
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class Cbcl2Text {

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) {
		//Parse arguments
		String runFolder=null;
		String outFile=null;
		int lane=-1;
		String tiles=null;
		String length=null;

		for(String arg : args){
			if(arg.indexOf('=')<0){continue;}
			String a=arg.split("=")[0].toLowerCase();
			String b=arg.substring(arg.indexOf('=')+1);

			if(a.equals("runfolder") || a.equals("in")){
				runFolder=b;
			} else if(a.equals("out") || a.equals("output")){
				outFile=b;
			} else if(a.equals("lane")){
				lane=Integer.parseInt(b);
			} else if(a.equals("tiles")){
				tiles=b;
			} else if(a.equals("length") || a.equals("lengths")){
				length=b;
			}
		}

		//Validate required arguments
		if(runFolder==null || outFile==null || lane<0){
			System.err.println("Usage: Cbcl2Text runfolder=<path> out=<file> lane=<int> [tiles=<list>] [length=auto|151,19,10,151]");
			System.err.println("Example: Cbcl2Text runfolder=./151T8B8B151T_cbcl out=output.txt lane=1");
			System.err.println("  length=auto - parse RunInfo.xml for read structure");
			System.err.println("  length=151,19,10,151 - manual read lengths");
			System.exit(1);
		}

		//Run conversion
		Cbcl2Text converter=new Cbcl2Text(runFolder, outFile, lane, tiles, length);
		converter.process();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	public Cbcl2Text(String runFolder_, String outFile_, int lane_, String tiles_, String length_) {
		runFolder=runFolder_;
		outFile=outFile_;
		lane=lane_;

		//Parse tile list if provided
		if(tiles_!=null && !tiles_.isEmpty()){
			String[] parts=tiles_.split(",");
			tileList=new IntList();
			for(String p : parts){
				tileList.add(Integer.parseInt(p.trim()));
			}
		}

		//Parse read lengths
		if(length_!=null && !length_.isEmpty()){
			if(length_.equals("auto")){
				//Parse RunInfo.xml
				try {
					readLengths=RunInfoParser.parseReadLengths(runFolder_);
					splitReads=true;
					System.err.println("Read structure from RunInfo.xml: " + Arrays.toString(readLengths));
				} catch(IOException e){
					System.err.println("Warning: Failed to parse RunInfo.xml: " + e.getMessage());
					System.err.println("Continuing with concatenated output");
				}
			}else{
				//Manual comma-separated lengths
				String[] parts=length_.split(",");
				readLengths=new int[parts.length];
				for(int i=0; i<parts.length; i++){
					readLengths[i]=Integer.parseInt(parts[i].trim());
				}
				splitReads=true;
				System.err.println("Read structure (manual): " + Arrays.toString(readLengths));
			}
		}

		//Build paths
		intensitiesDir=runFolder+"/Data/Intensities";
		baseCallsDir=intensitiesDir+"/BaseCalls";
		locsFile=intensitiesDir+"/s.locs";
	}

	/*--------------------------------------------------------------*/
	/*----------------       Core Processing        ----------------*/
	/*--------------------------------------------------------------*/

	public void process() {
		//Read cluster positions from s.locs
		System.err.println("Reading positions from " + locsFile);
		float[][] positions;
		try{
			positions=LocsReader.readPositions(locsFile);

			totalClusters=positions.length;
			System.err.println("Total clusters: " + totalClusters);

			//Determine which tiles to process
			determineTiles();

			System.err.println("\nWriting output to " + outFile);
			FileFormat ff=FileFormat.testOutput(outFile, FileFormat.FASTQ, null, true, true, false, false);
			ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
			long written=0;

			//Process each tile
			for(int i=0; i<tilesToProcess.size; i++){
				final int tileNum=tilesToProcess.get(i);
				System.err.println("\nProcessing tile " + tileNum);
				List<ClusterData> tileClusters=processTile(lane, tileNum, positions);
				written+=writeOutput(tileClusters, bsw, ff.fastq());
			}

			boolean error=bsw.poisonAndWait();
			if(error) {System.err.println("Something went wrong closing "+ff);}
			System.err.println("Complete! Wrote " + written + " clusters");
		}catch(IOException e){
			e.printStackTrace();
			KillSwitch.kill();
		}
	}

	/**
	 * Determine which tiles to process based on available files.
	 */
	private void determineTiles() throws IOException {
		tilesToProcess=new IntList();

		//If user specified tiles, use those
		if(tileList!=null && !tileList.isEmpty()){
			tilesToProcess.addAll(tileList);
			return;
		}

		//Otherwise, scan for available filter files
		String laneDir=baseCallsDir+"/L"+String.format("%03d", lane);
		File dir=new File(laneDir);
		File[] files=dir.listFiles();
		if(files==null){
			throw new IOException("Lane directory not found: " + laneDir);
		}

		//Look for s_1_TTTT.filter files
		for(File f : files){
			String name=f.getName();
			if(name.startsWith("s_1_") && name.endsWith(".filter")){
				//Extract tile number
				String tileStr=name.substring(4, name.length()-7);
				int tileNum=Integer.parseInt(tileStr);
				tilesToProcess.add(tileNum);
			}
		}

		tilesToProcess.sort();
		System.err.println("Found " + tilesToProcess.size() + " tiles: " + tilesToProcess);
	}

	/**
	 * Process a single tile across all cycles.
	 */
	private List<ClusterData> processTile(int lane, int tileNum, float[][] positions) throws IOException {
		//Read pass-filter flags
		String filterFile=baseCallsDir+"/L"+String.format("%03d", lane)+
		                  "/s_1_"+tileNum+".filter";
		boolean[] passFilter=FilterReader.readFilters(filterFile);
		int numClusters=passFilter.length;

		System.err.println("  Tile has " + numClusters + " clusters");

		//Determine read structure by scanning available cycles
		List<Integer> availableCycles=findAvailableCycles();
		if(availableCycles.isEmpty()){
			throw new IOException("No cycle directories found");
		}

		System.err.println("  Found " + availableCycles.size() + " cycles");

		//Create ClusterData objects
		List<ClusterData> clusters=new ArrayList<>();
		for(int i=0; i<numClusters; i++){
			float x=positions[i][0];
			float y=positions[i][1];
			ClusterData cd=new ClusterData(lane, tileNum, i, x, y);
			cd.passFilter=passFilter[i];
			clusters.add(cd);
		}

		//Read base calls from each cycle
		//For now, concatenate all cycles into a single read
		//TODO: Proper read structure parsing (151-8-8-151)
		ByteBuilder[] basesPerCluster=new ByteBuilder[numClusters];
		ByteBuilder[] qualsPerCluster=new ByteBuilder[numClusters];
		for(int i=0; i<numClusters; i++){
			basesPerCluster[i]=new ByteBuilder();
			qualsPerCluster[i]=new ByteBuilder();
		}

		//Determine which surface this tile is on by checking CBCL headers
		int surface=determineSurface(tileNum, availableCycles.get(0));
		System.err.println("  Tile " + tileNum + " is on surface " + surface);

		for(int cycle : availableCycles){
			String cbclFile=getCbclFilename(cycle, surface);
			if(!new File(cbclFile).exists()){
				System.err.println("  Warning: CBCL file not found for cycle " + cycle);
				continue;
			}

			//Decode this cycle
			byte[][] data=CbclDecoder.readTile(cbclFile, tileNum);
			byte[] bases=data[0];
			byte[] quals=data[1];

			//Append to each cluster's sequence
			for(int i=0; i<numClusters && i<bases.length; i++){
				basesPerCluster[i].append(bases[i]);
				qualsPerCluster[i].append(quals[i]);
			}
		}

		//Assign sequences to clusters - split by read if lengths specified
		for(int i=0; i<numClusters; i++){
			byte[] allBases=basesPerCluster[i].toBytes();
			byte[] allQuals=qualsPerCluster[i].toBytes();

			if(readLengths!=null && readLengths.length>0){
				//Split into reads based on specified lengths
				clusters.get(i).setData(allBases, allQuals, readLengths);
			}else{
				//Concatenated (default)
				clusters.get(i).basesR1=allBases;
				clusters.get(i).qualsR1=allQuals;
			}
		}

		return clusters;
	}

	/**
	 * Find all available cycle directories.
	 */
	private List<Integer> findAvailableCycles() {
		List<Integer> cycles=new ArrayList<>();
		String laneDir=baseCallsDir+"/L"+String.format("%03d", lane);
		File dir=new File(laneDir);
		File[] files=dir.listFiles();
		if(files==null){return cycles;}

		for(File f : files){
			String name=f.getName();
			if(name.startsWith("C") && name.endsWith(".1") && f.isDirectory()){
				//Extract cycle number
				String cycleStr=name.substring(1, name.length()-2);
				try {
					int cycle=Integer.parseInt(cycleStr);
					cycles.add(cycle);
				} catch(NumberFormatException e){
					//Skip non-numeric directories
				}
			}
		}

		Collections.sort(cycles);
		return cycles;
	}

	/**
	 * Determine which surface (1 or 2) contains the specified tile.
	 */
	private int determineSurface(int tileNum, int sampleCycle) throws IOException {
		//Check surface 1
		String cbclFile1=getCbclFilename(sampleCycle, 1);
		if(new File(cbclFile1).exists()){
			CbclHeader header=new CbclHeader(cbclFile1);
			if(header.tileMetadata.containsKey(tileNum)){
				return 1;
			}
		}

		//Check surface 2
		String cbclFile2=getCbclFilename(sampleCycle, 2);
		if(new File(cbclFile2).exists()){
			CbclHeader header=new CbclHeader(cbclFile2);
			if(header.tileMetadata.containsKey(tileNum)){
				return 2;
			}
		}

		throw new IOException("Tile " + tileNum + " not found in any surface CBCL file");
	}

	/**
	 * Get CBCL filename for a specific cycle and surface.
	 */
	private String getCbclFilename(int cycle, int surface) {
		//Note: directory name has .1 suffix (always), but inside are surface 1 and 2 files
		return baseCallsDir+"/L"+String.format("%03d", lane)+
		       "/C"+cycle+".1/L"+String.format("%03d", lane)+"_"+surface+".cbcl";
	}

	/**
	 * Write output to tab-delimited text file or fastq.
	 */
	private long writeOutput(List<ClusterData> clusters, ByteStreamWriter bsw, boolean fastq) {
		ByteBuilder bb=new ByteBuilder();
		long written=0;
		for(ClusterData cd : clusters){
			if(cd.passFilter || !passOnly) {
				bb.clear();
				if(fastq) {cd.toBytes(bb);}
				else{cd.toBytesOld(bb, splitReads);}
				bsw.print(bb);
				written++;
			}
		}
		return written;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	private final String runFolder;
	private String outFile="stdout.fastq";
	private final int lane;
	private IntList tileList=null;

	private String intensitiesDir;
	private String baseCallsDir;
	private String locsFile;

	private int totalClusters;
	private IntList tilesToProcess;

	private int[] readLengths=null; //Read lengths for splitting
	private boolean splitReads=false; //Output comma-delimited reads
	private boolean passOnly=true; //Output comma-delimited reads
}
