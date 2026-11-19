package prot.cluster;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;

/** Partitions all protein input ordinals into near-equal seeded blocks, without clustering.
 * Uses the native argument and preserved FASTA loading paths of DdlClust.
 * The optional hydrophobicity ordering is stable after shuffling, so ties are seeded.
 * in=file1,file2 preserves declared file order; tailcrop=N crops before sorting.
 * ready=part#.ready optionally signals each closed shard, before whole-map completion.
 * @author Keqing
 */
public final class DdlClustPartition{

	public static void main(String[] args){
		final DdlClustPartition tool=new DdlClustPartition(args);
		try{tool.process();}finally{Shared.closeStream(tool.outstream);}
	}

	/** Validates source/output identities before any writer can truncate the source. */
	public DdlClustPartition(String[] args){
		final PreParser pre=new PreParser(args,getClass(),false);
		args=pre.args; outstream=pre.outstream;
		final Parser parser=new Parser(); parser.overwrite=false;
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] pair=Parse.splitOnFirst(arg,'=');
			final String key=pair[0].toLowerCase(Locale.ROOT);
			final String value=pair.length>1 && !pair[1].equalsIgnoreCase("null") ? pair[1] : null;
			if(key.equals("ways")){ways=Parse.parseIntKMG(value);}
			else if(key.equals("tailcrop")){tailCrop=Parse.parseIntKMG(value);}
			else if(key.equals("map")){mapping=value;}
			else if(key.equals("ready")){readyPattern=value;}
			else if(key.equals("seed")){seed=Long.parseLong(value);}
			else if(key.equals("order")){
				if(value==null){throw new IllegalArgumentException("order requires shuffle or hydro");}
				order=ClusteringOrder.Mode.valueOf(value.toUpperCase(Locale.ROOT));
			}else if(parser.parse(arg,key,value)){}
			else if(i==0 && Tools.looksLikeInputSequenceStream(arg)){parser.in1=arg;}
			else{throw new IllegalArgumentException("Unknown partition parameter: "+arg);}
		}
		output=parser.out1; overwrite=parser.overwrite;
		if(parser.in1==null || output==null || mapping==null || parser.in2!=null || parser.out2!=null ||
			parser.append || parser.maxReads>=0 || ways<1 || tailCrop<0 ||
			(order!=ClusteringOrder.Mode.SHUFFLE && order!=ClusteringOrder.Mode.HYDRO)){
			throw new IllegalArgumentException("Partitioning requires complete FASTA inputs, out=, map=, ways>0, tailcrop>=0 and order=shuffle|hydro; paired input, append and reads limits are unsupported");
		}
		if(output.indexOf('#')<0 || output.indexOf('#')!=output.lastIndexOf('#')){
			throw new IllegalArgumentException("Partition out= requires exactly one # placeholder");
		}
		if(readyPattern!=null && (readyPattern.indexOf('#')<0 || readyPattern.indexOf('#')!=readyPattern.lastIndexOf('#'))){
			throw new IllegalArgumentException("Partition ready= requires exactly one # placeholder");
		}
		inputs=ClusteringInput.inputFiles(parser.in1);
		for(String input:inputs){
			final Path source=Paths.get(input);
			if(!Files.isRegularFile(source) || !Files.isReadable(source)){
				throw new IllegalArgumentException("Partition input must be a readable disk FASTA: "+input);
			}
			checkSourceAlias(source,mapping);
			for(int i=0; i<ways; i++){
				checkSourceAlias(source,output.replace("#",Integer.toString(i)));
				if(readyPattern!=null){checkSourceAlias(source,readyPattern.replace("#",Integer.toString(i)));}
			}
		}
	}

	/** Native loading retains every ordinal, including short proteins and duplicate text IDs. */
	private void process(){
		Shared.AMINO_IN=true; Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=Read.LOWER_CASE_TO_N=Read.DOT_DASH_X_TO_N=Read.FIX_HEADER=false;
		Read.JUNK_MODE=Read.IGNORE_JUNK;
		final ClusteringInput loaded=ClusteringInput.load(inputs,-1,true,tailCrop);
		ClusteringShardWriter.write(loaded,output,mapping,ways,order,seed,overwrite,readyPattern);
		final int minimum=loaded.sequences.length/ways, extra=loaded.sequences.length%ways;
		assert((long)minimum*ways+extra==loaded.sequences.length) : "Balanced shard counts must conserve every global input ordinal";
		outstream.println("Partitioned "+loaded.sequences.length+" proteins into "+ways+" files of "+
			minimum+" records plus one in the first "+extra+" files; order="+
			order.toString().toLowerCase(Locale.ROOT)+", seed="+seed+", tailcrop="+tailCrop+".");
	}

	/** Checking only source-to-output pairs is linear in requested files, even for hard links. */
	private static void checkSourceAlias(Path source,String destination){
		try{
			final Path path=Paths.get(destination);
			if(Files.exists(path) && Files.isSameFile(source,path)){
				throw new IllegalArgumentException("Partition output aliases the source FASTA: "+destination);
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not check partition input/output identity",failure);}
	}

	private final PrintStream outstream;
	private final String[] inputs;
	private final String output;
	private final boolean overwrite;
	private String mapping;
	private String readyPattern;
	private int ways=4;
	private int tailCrop=0;
	private long seed=1;
	private ClusteringOrder.Mode order=ClusteringOrder.Mode.SHUFFLE;
}
