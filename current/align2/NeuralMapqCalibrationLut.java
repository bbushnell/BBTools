package align2;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

/** Immutable monotone raw-error to calibrated-error lookup table. */
public final class NeuralMapqCalibrationLut {

	public static NeuralMapqCalibrationLut load(final String path,final int mapqCap)throws IOException{
		if(mapqCap<0||mapqCap>60){throw new IllegalArgumentException("MAPQ cap outside [0,60]: "+mapqCap);}
		final ArrayList<Float> highs=new ArrayList<Float>();final ArrayList<Double> probabilities=new ArrayList<Double>();
		try(BufferedReader reader=Files.newBufferedReader(Paths.get(path),StandardCharsets.UTF_8)){
			final String header=reader.readLine();
			if(!HEADER.equals(header)){throw new IllegalArgumentException("Unexpected neural MAPQ LUT header in "+path+": "+header);}
			String line;float previousHigh=Float.NEGATIVE_INFINITY;double previousProbability=-1;
			while((line=reader.readLine())!=null){
				if(line.isEmpty()){continue;}final String[] fields=line.split("\t",-1);
				if(fields.length!=7){throw new IllegalArgumentException("Neural MAPQ LUT row width differs: "+line);}
				final float low=Float.parseFloat(fields[1]),high=Float.parseFloat(fields[2]);
				final long rows=Long.parseLong(fields[3]),errors=Long.parseLong(fields[4]);
				final double probability=Double.parseDouble(fields[5]);
				if(!Float.isFinite(low)||!Float.isFinite(high)||low>high||low<previousHigh||rows<1||errors<0||errors>rows||
						!Double.isFinite(probability)||probability<=0||probability>1||probability<previousProbability){
					throw new IllegalArgumentException("Invalid or nonmonotone neural MAPQ LUT row: "+line);
				}
				highs.add(high);probabilities.add(probability);previousHigh=high;previousProbability=probability;
			}
		}
		if(highs.isEmpty()){throw new IllegalArgumentException("Neural MAPQ LUT has no rows: "+path);}
		final float[] h=new float[highs.size()];final double[] p=new double[probabilities.size()];
		for(int i=0;i<h.length;i++){h[i]=highs.get(i);p[i]=probabilities.get(i);}
		return new NeuralMapqCalibrationLut(h,p,mapqCap,path);
	}

	private NeuralMapqCalibrationLut(final float[] high_,final double[] probability_,final int cap_,final String source_){
		high=high_;probability=probability_;mapqCap=cap_;source=source_;
	}

	public double calibratedError(final float rawError){
		if(!Float.isFinite(rawError)||rawError<0||rawError>1){throw new IllegalArgumentException("Raw neural error outside [0,1]: "+rawError);}
		int low=0,hi=high.length-1;while(low<hi){final int mid=(low+hi)>>>1;if(rawError<=high[mid]){hi=mid;}else{low=mid+1;}}
		return probability[low];
	}

	public int mapq(final float rawError){
		return mapq(rawError,mapqCap);
	}
	public int mapq(final float rawError,final int cap){
		if(cap<0||cap>60){throw new IllegalArgumentException("MAPQ cap outside [0,60]: "+cap);}
		final double p=calibratedError(rawError);final int q=(int)Math.round(-10*Math.log10(p));return Math.min(cap,Math.max(0,q));
	}

	public int size(){return high.length;}
	public int mapqCap(){return mapqCap;}
	public String source(){return source;}

	private final float[] high;private final double[] probability;private final int mapqCap;private final String source;
	private static final String HEADER="block\traw_low\traw_high\trows\terrors\tcalibrated_error\tcalibrated_mapq";
}
