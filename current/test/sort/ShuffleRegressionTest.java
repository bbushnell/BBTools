package test.sort;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.Parse;

/** Small deterministic FASTQ fixtures and independent record/mate checks for the shuffle tools.
 * @author Shinobu
 */
public class ShuffleRegressionTest{

	/** Runs generate=directory or check=file; format=fa checks one-line FASTA and order=sequence checks fixture ordering. */
	public static void main(String[] args){
		String generate=null, check=null, check2=null;
		int count=512, mate=0;
		boolean fasta=false, sequenceOrder=false;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			if(split[0].equals("generate")){generate=split[1];}
			else if(split[0].equals("check")){check=split[1];}
			else if(split[0].equals("check2")){check2=split[1];}
			else if(split[0].equals("count")){count=Integer.parseInt(split[1]);}
			else if(split[0].equals("mate")){mate=Integer.parseInt(split[1]);}
			else if(split[0].equals("order")){
				if(!split[1].equals("sequence")){throw new IllegalArgumentException("Expected order=sequence");}
				sequenceOrder=true;
			}else if(split[0].equals("format")){
				if(!split[1].equals("fa") && !split[1].equals("fq")){
					throw new IllegalArgumentException("Expected format=fa or fq");
				}
				fasta=split[1].equals("fa");
			}else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(count<0 || mate<0 || mate>2 || (generate==null)==(check==null)){
			throw new IllegalArgumentException("Choose generate or check; count>=0 and mate=0..2");
		}
		if(generate!=null){generate(generate, count);}
		else{
			final int[] order=check(check, count, mate, fasta);
			if(sequenceOrder){checkSequenceOrder(order, mate);}
			if(check2!=null){
				if(mate!=1){throw new IllegalArgumentException("check2 requires mate=1");}
				final int[] mateOrder=check(check2, count, 2, fasta);
				if(!Arrays.equals(order, mateOrder)){throw new IllegalStateException("Mate output orders differ");}
				System.out.println("PASS\tpaired-order\t"+count);
			}
		}
	}

	/** Checks this fixture's equal-length, unique sequences with Java's independent String ordering. */
	private static void checkSequenceOrder(int[] order, int mate){
		String previous=null;
		for(int id : order){
			final String current=new String(sequence(id, mate), StandardCharsets.US_ASCII);
			if(previous!=null && previous.compareTo(current)>0){
				throw new IllegalStateException("Sequence order decreased: "+previous+" > "+current);
			}
			previous=current;
		}
		System.out.println("PASS\tsequence-order\t"+order.length);
	}

	/** Creates three known-order inputs; each read's sequence identifies its numeric ID and mate. */
	private static void generate(String directory, int count){
		final File dir=new File(directory);
		if(!dir.isDirectory() && !dir.mkdirs()){throw new IllegalStateException("Cannot create "+dir);}
		for(int mate=0; mate<=2; mate++){
			final String name=mate==0 ? "single.fq" : "r"+mate+".fq";
			final ByteStreamWriter out=new ByteStreamWriter(new File(dir, name).toString(), false, false, false);
			out.start();
			for(int id=0; id<count; id++){
				out.print("@read").print(id);
				if(mate>0){out.print('/').print(mate);}
				out.nl().println(sequence(id, mate)).println("+").println(QUALITY);
			}
			if(out.poisonAndWait()){throw new IllegalStateException("Fixture write failed: "+name);}
		}
	}

	/** Checks exact multiplicity, labels, bases and FASTQ qualities, returning record ID order. */
	private static int[] check(String path, int count, int mate, boolean fasta){
		if(!new File(path).isFile()){throw new IllegalStateException("Missing output: "+path);}
		final boolean[] seen=new boolean[count];
		final int[] order=new int[count];
		final ByteFile in=ByteFile.makeByteFile(path, false);
		int records=0;
		try{
			for(byte[] header=in.nextLine(); header!=null; header=in.nextLine()){
				final String name=new String(header, StandardCharsets.US_ASCII);
				if(!name.startsWith(fasta ? ">read" : "@read")){throw new IllegalStateException("Unexpected ID: "+name);}
				final int slash=name.indexOf('/');
				if(mate==0 ? slash>=0 : !name.endsWith("/"+mate)){
					throw new IllegalStateException("Wrong mate in "+path+": "+name);
				}
				final int id=Parse.parseInt(name, 5, slash<0 ? name.length() : slash);
				if(id<0 || id>=count || seen[id]){throw new IllegalStateException("Duplicate/invalid ID: "+name);}
				seen[id]=true;
				final byte[] bases=in.nextLine();
				if(!Arrays.equals(bases, sequence(id, mate))){
					throw new IllegalStateException("Changed or truncated sequence: "+name);
				}
				if(!fasta && (!Arrays.equals(in.nextLine(), PLUS) || !Arrays.equals(in.nextLine(), QUALITY))){
					throw new IllegalStateException("Changed or truncated FASTQ quality: "+name);
				}
				order[records++]=id;
			}
		}finally{
			if(in.close()){throw new IllegalStateException("Input failure: "+path);}
		}
		if(records!=count){throw new IllegalStateException("Expected "+count+" records; found "+records);}
		System.out.println("PASS\t"+path+"\t"+records+"\tmate="+mate);
		return order;
	}

	/** Encodes ID in the first bases and mate in the last base, independently of output order. */
	private static byte[] sequence(int id, int mate){
		assert(id>=0 && mate>=0 && mate<=2) : "Fixture IDs/mates must fit their independent oracle";
		final byte[] bases=new byte[20];
		for(int i=0; i<bases.length; i++){bases[i]=ALPHABET[(id>>>(2*(i&15)))&3];}
		bases[bases.length-1]=ALPHABET[mate];
		return bases;
	}

	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
	private static final byte[] PLUS={'+'};
	private static final byte[] QUALITY="IIIIIIIIIIIIIIIIIIII".getBytes(StandardCharsets.US_ASCII);
}
