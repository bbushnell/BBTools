package align2;

import fileIO.TextFile;

/**
 * Legacy partial extractor for Mapping timing lines and GradeSamFile report names.
 * Emits mapping time in seconds followed by the experiment name, under a
 * matching two-column header. It does not extract accuracy metrics.
 * Requires each Mapping timing line to precede its matching report name.
 * A row is written only after both fields are validated; missing names/times and
 * duplicate pending timing lines are rejected. Elapsed lines are ignored.
 * Reads the entire input into memory and writes to stdout.
 *
 * @author Brian Bushnell
 * @date 2014
 */
public class ReformatBatchOutput2{

	/**
	 * Program entry point that processes a mapping statistics file.
	 * Reads the input file specified as the first command-line argument,
	 * emits mapping-time/name rows for complete historical timing/name pairs.
	 * @param args Command-line arguments where args[0] is the input file path
	 */
	public static void main(String[] args){
		assert(args.length>0) : "A legacy mapping-report input path is required";
		TextFile tf=new TextFile(args[0], false);
		final String[] lines;
		try{lines=tf.toStringLines();}
		finally{tf.close();}

		String mappingTime=null;

		System.out.println(header());

		for(String s : lines){
			if(s.startsWith("Mapping Statistics for ")){
				final String name=s.substring("Mapping Statistics for ".length()).replace(".sam:", "").trim();
				if(mappingTime==null || name.isEmpty() || name.indexOf('\t')>=0){
					throw new IllegalArgumentException("A report name requires a preceding Mapping time and one nonempty name field: "+s);
				}
				System.out.println(mappingTime+"\t"+name);
				mappingTime=null;
			}else if(s.startsWith("Mapping:")){
				if(mappingTime!=null){throw new IllegalArgumentException("Mapping time has no report name before the next timing line: "+s);}
				final String value=s.substring("Mapping:".length()).replace("seconds.", "").trim();
				final double seconds;
				try{seconds=Double.parseDouble(value);}
				catch(NumberFormatException e){throw new IllegalArgumentException("Invalid Mapping time: "+s, e);}
				if(seconds<0 || Double.isNaN(seconds) || Double.isInfinite(seconds)){
					throw new IllegalArgumentException("Mapping time must be finite and nonnegative: "+s);
				}
				mappingTime=value;
			}
		}
		if(mappingTime!=null){throw new IllegalArgumentException("Mapping time has no report name at end of input: "+mappingTime);}
	}

	/**
	 * Returns the tab-separated header string for output formatting.
	 * Names the emitted mapping time (seconds) and experiment name, in that order.
	 * @return Tab-separated header string with column names
	 */
	public static String header(){
		return "mapTime\tname";
	}

}
