package clade;

import java.util.ArrayList;

import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetParser;

/**
 * Loads neural networks and calibration constants from a .bbnets file.
 * Supports hybrid format: ##network L B where B=-1 means all-length,
 * B>=0 means per-length-bin.
 *
 * @author Ady
 */
public class SerialNNLoader {

	public static LoadedNets load(String path) {
		if(path==null){return null;}
		ArrayList<byte[]> allLines=ByteFile.toLines(path);
		if(allLines==null || allLines.isEmpty()){return null;}

		int levels=-1, bins=-1;

		int pos=0;
		while(pos<allLines.size()) {
			String line=new String(allLines.get(pos));
			if(line.startsWith("#levels ")){
				levels=Integer.parseInt(line.substring(8).trim());
			}else if(line.startsWith("#bins ")){
				bins=Integer.parseInt(line.substring(6).trim());
			}else if(line.startsWith("##network ")){
				break;
			}
			pos++;
		}

		if(levels<1){return null;}
		if(bins<1){bins=1;}
		CellNet[] allLenNets=new CellNet[levels];
		float[][] allLenCal=new float[levels][4];
		float[][] allLenLutX=new float[levels][], allLenLutY=new float[levels][];
		String[] allLenLabels=new String[levels];
		CellNet[][] binNets=new CellNet[levels][bins];
		float[][][] binCal=new float[levels][bins][4];
		float[][][] binLutX=new float[levels][bins][], binLutY=new float[levels][bins][];

		while(pos<allLines.size()) {
			String line=new String(allLines.get(pos));
			if(!line.startsWith("##network ")){pos++; continue;}

			String[] parts=line.substring(10).trim().split("\\s+");
			int lvl=Integer.parseInt(parts[0]);
			int bin=parts.length>1 ? Integer.parseInt(parts[1]) : -1;
			pos++;

			float[] calParams=null, lutX=null, lutY=null;
			String label=null;
			ArrayList<byte[]> netLines=new ArrayList<>();

			while(pos<allLines.size()) {
				String s=new String(allLines.get(pos));
				if(s.startsWith("##network ")){break;}
				if(s.startsWith("#cal ")){
					String[] cp=s.substring(5).trim().split("\\s+");
					calParams=new float[]{
						Float.parseFloat(cp[0]), Float.parseFloat(cp[1]),
						Float.parseFloat(cp[2]), Float.parseFloat(cp[3])
					};
				}else if(s.startsWith("#lut ")){
					//Isotonic calibration table: a count, then exactly that many "x<tab>y" lines.
					//Consumed BY COUNT because the points are bare numbers that CellNetParser would
					//otherwise ingest as network lines.
					final int n=Integer.parseInt(s.substring(5).trim());
					assert(n>0 && pos+n<allLines.size()) : "#lut "+n+" overruns the file at line "+pos;
					lutX=new float[n]; lutY=new float[n];
					for(int i=0; i<n; i++){
						pos++;
						final String[] xy=new String(allLines.get(pos)).trim().split("\\s+");
						assert(xy.length>=2) : "Malformed #lut point at line "+pos;
						lutX[i]=Float.parseFloat(xy[0]); lutY[i]=Float.parseFloat(xy[1]);
						assert(i==0 || lutX[i]>=lutX[i-1]) : "#lut x must be ascending at line "+pos;
					}
				}else if(s.startsWith("#label")){
					//The bundle DECLARES which taxonomic level each network answers, so the
					//consumer never has to assume a fixed level list. Dropping a level from the
					//file is then sufficient to remove it everywhere.
					label=s.substring(6).trim();
				}else{
					netLines.add(allLines.get(pos));
				}
				pos++;
			}

			if(lvl<0 || lvl>=levels){continue;}
			CellNet net=netLines.isEmpty() ? null : CellNetParser.loadFromLines(netLines);
			//Trusted internal-file assumptions (.bbnets is generated, not user input): a bin>=bins is silently dropped (neither branch stores it); a net with no #cal keeps the all-zero default cal -> CladeConfidence.calibrate returns 0 (effectively no-confidence). Both out-of-scope per input-validity (malformed generated file), but noted. A null net (empty block) is handled gracefully downstream (sigmoid fallback).
			if(bin<0){
				allLenNets[lvl]=net;
				if(calParams!=null){allLenCal[lvl]=calParams;}
				if(lutX!=null){allLenLutX[lvl]=lutX; allLenLutY[lvl]=lutY;}
				allLenLabels[lvl]=label;
			}else if(bin<bins){
				binNets[lvl][bin]=net;
				if(calParams!=null){binCal[lvl][bin]=calParams;}
				if(lutX!=null){binLutX[lvl][bin]=lutX; binLutY[lvl][bin]=lutY;}
			}
		}

		LoadedNets result=new LoadedNets();
		result.allLenNets=allLenNets;
		result.allLenCal=allLenCal;
		result.allLenLutX=allLenLutX;
		result.allLenLutY=allLenLutY;
		result.allLenLabels=allLenLabels;
		result.binNets=binNets;
		result.binCal=binCal;
		result.binLutX=binLutX;
		result.binLutY=binLutY;
		result.levels=levels;
		result.bins=bins;
		return result;
	}

	public static class LoadedNets {
		public CellNet[] allLenNets;
		public float[][] allLenCal;
		/** Isotonic calibration table per level (null when the bundle carries only #cal). */
		public float[][] allLenLutX, allLenLutY;
		/** Level name each network was trained for, from its #label line; null if absent. */
		public String[] allLenLabels;
		public CellNet[][] binNets;
		public float[][][] binCal;
		public float[][][] binLutX, binLutY;
		public int levels;
		public int bins;
	}

}
