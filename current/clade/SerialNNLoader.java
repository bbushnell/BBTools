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
		CellNet[][] binNets=new CellNet[levels][bins];
		float[][][] binCal=new float[levels][bins][4];

		while(pos<allLines.size()) {
			String line=new String(allLines.get(pos));
			if(!line.startsWith("##network ")){pos++; continue;}

			String[] parts=line.substring(10).trim().split("\\s+");
			int lvl=Integer.parseInt(parts[0]);
			int bin=parts.length>1 ? Integer.parseInt(parts[1]) : -1;
			pos++;

			float[] calParams=null;
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
				}else if(s.startsWith("#label")){
					// skip
				}else{
					netLines.add(allLines.get(pos));
				}
				pos++;
			}

			if(lvl<0 || lvl>=levels){continue;}
			CellNet net=netLines.isEmpty() ? null : CellNetParser.loadFromLines(netLines);
			if(bin<0){
				allLenNets[lvl]=net;
				if(calParams!=null){allLenCal[lvl]=calParams;}
			}else if(bin<bins){
				binNets[lvl][bin]=net;
				if(calParams!=null){binCal[lvl][bin]=calParams;}
			}
		}

		LoadedNets result=new LoadedNets();
		result.allLenNets=allLenNets;
		result.allLenCal=allLenCal;
		result.binNets=binNets;
		result.binCal=binCal;
		result.levels=levels;
		result.bins=bins;
		return result;
	}

	public static class LoadedNets {
		public CellNet[] allLenNets;
		public float[][] allLenCal;
		public CellNet[][] binNets;
		public float[][][] binCal;
		public int levels;
		public int bins;
	}

}
