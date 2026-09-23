package align2;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

/** Immutable V2 MAPQ caps indexed by mode, reference regime, and read length. */
public final class NeuralMapqLengthCaps {

	public static NeuralMapqLengthCaps load(final String path)throws IOException{
		final int[][][] caps=new int[2][2][BANDS.length];for(int m=0;m<caps.length;m++)for(int r=0;r<caps[m].length;r++)for(int b=0;b<caps[m][r].length;b++){caps[m][r][b]=-1;}
		int rows=0;
		try(BufferedReader reader=Files.newBufferedReader(Paths.get(path),StandardCharsets.UTF_8)){
			final String header=reader.readLine();if(!HEADER.equals(header)){throw new IllegalArgumentException("Unexpected neural MAPQ cap header in "+path+": "+header);}
			String line;while((line=reader.readLine())!=null){
				if(line.isEmpty()){continue;}final String[] f=line.split("\t",-1);if(f.length!=9){throw new IllegalArgumentException("Neural MAPQ cap row width differs: "+line);}
				final int mode=modeIndex(f[0]),regime=regimeIndex(f[1]),band=bandIndex(f[2]),cap=Integer.parseInt(f[3]),threshold=Integer.parseInt(f[5]);
				final long retained=Long.parseLong(f[6]),errors=Long.parseLong(f[7]);final double evidenceQ=Double.parseDouble(f[8]);
				if(mode<0||regime<0||band<0||cap<0||cap>50||f[4].isEmpty()||threshold<0||threshold>50||retained<1||errors<0||errors>retained||!Double.isFinite(evidenceQ)||cap!=(int)Math.floor(evidenceQ)||caps[mode][regime][band]>=0){throw new IllegalArgumentException("Invalid or duplicate neural MAPQ cap row: "+line);}
				caps[mode][regime][band]=cap;rows++;
			}
		}
		if(rows!=28){throw new IllegalArgumentException("Expected28 neural MAPQ caps, observed "+rows);}
		for(int m=0;m<caps.length;m++)for(int r=0;r<caps[m].length;r++)for(int b=0;b<caps[m][r].length;b++)if(caps[m][r][b]<0){throw new IllegalArgumentException("Missing neural MAPQ cap mode="+m+" regime="+r+" band="+b);}
		return new NeuralMapqLengthCaps(caps,path);
	}

	private NeuralMapqLengthCaps(final int[][][] caps_,final String source_){caps=caps_;source=source_;}

	/** Returns -1 when V2 must fall back to legacy MAPQ. */
	public int cap(final boolean paired,final int referenceRegime,final int readLength){
		final int regime=referenceRegime==NeuralMapqReferenceScale.LARGE?0:referenceRegime==NeuralMapqReferenceScale.SMALL?1:-1;
		final int band=bandForLength(readLength);if(regime<0||band<0){return -1;}return caps[paired?1:0][regime][band];
	}

	public String source(){return source;}
	public static int bandForLength(final int length){if(length<MINIMUM_LENGTH||length>MAXIMUM_LENGTH){return -1;}for(int i=0;i<BAND_ENDS.length;i++)if(length<=BAND_ENDS[i]){return i;}throw new AssertionError("Supported V2 length must match a band: "+length);}
	private static int modeIndex(final String value){return value.equals("single")?0:value.equals("paired")?1:-1;}
	private static int regimeIndex(final String value){return value.equals("large")?0:value.equals("small")?1:-1;}
	private static int bandIndex(final String value){for(int i=0;i<BANDS.length;i++)if(BANDS[i].equals(value)){return i;}return -1;}

	private final int[][][] caps;private final String source;
	public static final int MINIMUM_LENGTH=50,MAXIMUM_LENGTH=250;
	private static final int[] BAND_ENDS={74,99,124,174,199,224,250};
	private static final String[] BANDS={"50-74","75-99","100-124","125-174","175-199","200-224","225-250"};
	private static final String HEADER="mode\tregime\tband\tinteger_cap\tbinding_domain\tnominal_threshold\tretained\terrors\tevidence_mapq";
}
