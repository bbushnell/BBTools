package prok;

import java.nio.charset.StandardCharsets;

import fileIO.ByteFile;
import map.ObjectIntMap;

/** Per-model HBM endpoint coverage constants for ncRNA boundary feature layout v2. */
public final class NcrnaBoundaryFuzzConstants {

	static NcrnaBoundaryFuzzConstants load(String path, String[] modelNames){
		if(path==null){throw new IllegalArgumentException("boundaryfeatures=v2 requires boundaryfuzz=<model/end constants TSV>");}
		if(modelNames==null || modelNames.length<1){throw new IllegalArgumentException("Cannot bind boundary fuzz constants without model names");}
		final ObjectIntMap<String> index=new ObjectIntMap<String>(modelNames.length*2, String.class);
		for(int i=0; i<modelNames.length; i++){
			final String name=modelNames[i];
			if(name==null || name.length()<1){throw new IllegalArgumentException("Empty boundary model name at index "+i);}
			if(index.put(name,i)>=0){throw new IllegalArgumentException("Duplicate boundary model name: "+name);}
		}
		final float[][] start=new float[modelNames.length][3], stop=new float[modelNames.length][3];
		final boolean[][] cliff=new boolean[modelNames.length][2];
		final boolean[][] seen=new boolean[modelNames.length][2];
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length<1 || line[0]=='#'){continue;}
				final String s=new String(line,StandardCharsets.UTF_8);
				final String[] split=s.split("\\t",-1);
				if(split.length>1 && split[0].equalsIgnoreCase("model") && split[1].equalsIgnoreCase("end")){
					if(split.length==6 && split[2].equals("fuzz_minus") && split[3].equals("fuzz_center") && split[4].equals("fuzz_plus") && split[5].equals("cliff")){continue;}
					throw new IllegalArgumentException("Boundary fuzz header must be: model, end, fuzz_minus, fuzz_center, fuzz_plus, cliff: "+s);
				}
				if(split.length!=6){throw new IllegalArgumentException("Boundary fuzz row must have 6 tab-delimited fields: "+s);}
				final int model=index.get(split[0]);
				if(model<0){throw new IllegalArgumentException("Unknown boundary fuzz model '"+split[0]+"' in "+path);}
				final int end=parseEnd(split[1]);
				if(seen[model][end]){throw new IllegalArgumentException("Duplicate boundary fuzz row for "+split[0]+" end="+split[1]);}
				final float[] dest=(end==0 ? start[model] : stop[model]);
				for(int i=0; i<3; i++){
					dest[i]=Float.parseFloat(split[i+2]);
					if(!Float.isFinite(dest[i]) || dest[i]<0f || dest[i]>1f){throw new IllegalArgumentException(
						"Boundary fuzz values must be finite normalized coverage in [0,1]: "+s);}
				}
				cliff[model][end]=parseCliff(split[5],s);
				seen[model][end]=true;
			}
		}finally{
			if(bf.close()){throw new RuntimeException("Read error while loading boundary fuzz constants: "+path);}
		}
		for(int model=0; model<modelNames.length; model++){
			if(!seen[model][0] || !seen[model][1]){throw new IllegalArgumentException(
				"Boundary fuzz constants require start and stop rows for model "+modelNames[model]);}
		}
		return new NcrnaBoundaryFuzzConstants(start,stop,cliff);
	}

	float[] values(int model, boolean useStart){
		if(model<0 || model>=start.length){throw new IllegalArgumentException("Boundary fuzz model index out of range: "+model);}
		return useStart ? start[model] : stop[model];
	}

	boolean cliff(int model, boolean useStart){
		if(model<0 || model>=start.length){throw new IllegalArgumentException("Boundary fuzz model index out of range: "+model);}
		return cliff[model][useStart ? 0 : 1];
	}

	private NcrnaBoundaryFuzzConstants(float[][] start_,float[][] stop_,boolean[][] cliff_){start=start_;stop=stop_;cliff=cliff_;}

	private static int parseEnd(String s){
		if(s.equalsIgnoreCase("start") || s.equals("5")){return 0;}
		if(s.equalsIgnoreCase("stop") || s.equals("3")){return 1;}
		throw new IllegalArgumentException("Boundary fuzz end must be start/stop (or 5/3): "+s);
	}

	private static boolean parseCliff(String value,String row){
		if(value.equalsIgnoreCase("t") || value.equalsIgnoreCase("true")){return true;}
		if(value.equalsIgnoreCase("f") || value.equalsIgnoreCase("false")){return false;}
		throw new IllegalArgumentException("Boundary fuzz cliff must be t/f or true/false: "+row);
	}

	private final float[][] start,stop;
	private final boolean[][] cliff;
}
