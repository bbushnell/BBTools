package ml;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Utilities for manipulating neural-network training-vector files (the tab-delimited '#dims'-header
 * format read by ml.DataLoader / ml.Trainer). Combines any number of input vector files and emits one
 * or more outputs, applying — in this order — deduplicate, subsample, shuffle, partition (by fraction),
 * and positive/negative balance. The '#dims' header is carried through correctly (never shuffled into
 * the body), which is why this exists: a plain `shuf` corrupts the header.
 *
 * Example:
 *   vectorutils.sh a.tsv b.tsv c.tsv out=training.tsv:0.9,validation.tsv:0.1 shuffle samplerate=0.5 balance=0.3
 *
 * Params:
 *   in=a,b,c            Input vector files (also accepted as bare positional filenames).
 *   out=f1:0.9,f2:0.1   Output files, each with a partition fraction (fractions normalized to sum 1;
 *                       a single out= with no fraction gets the whole set).
 *   shuffle             Randomly shuffle rows (needed for a meaningful random partition).
 *   samplerate=X        Keep a random fraction X of rows (subsample). Default 1.
 *   balance=X           Upsample the MINORITY class (label = last column >= 0.5) with random duplicates
 *                       until it is fraction X of that output's total (Brian's default; X~0.25-0.3).
 *                       Applied to EACH output independently, after partition. Default 0 (off).
 *   deduplicate (dedupe) Remove exact-duplicate rows (sort-based). Default off.
 *   seed=N              RNG seed (>=0 deterministic; -1 random). Default 1.
 *   ow=t                Overwrite outputs.
 *
 * @author Noire, Brian Bushnell
 */
public class VectorUtils {

	public static void main(String[] args){
		Timer t=new Timer();
		VectorUtils x=new VectorUtils(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public VectorUtils(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("in") || a.equals("input")){
				for(String s : b.split(",")){if(s.length()>0){inFiles.add(s);}}
			}else if(a.equals("out") || a.equals("output")){
				parseOut(b);
			}else if(a.equals("shuffle")){
				shuffle=(b==null || Parse.parseBoolean(b));
			}else if(a.equals("samplerate") || a.equals("sample") || a.equals("subsample")){
				samplerate=Float.parseFloat(b);
			}else if(a.equals("balance")){
				balance=Float.parseFloat(b);
			}else if(a.equals("deduplicate") || a.equals("dedupe")){
				dedupe=(b==null || Parse.parseBoolean(b));
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("ow") || a.equals("overwrite")){
				overwrite=Parse.parseBoolean(b);
			}else if(arg.indexOf('=')<0){
				inFiles.add(arg);//bare positional filename
			}else{
				outstream.println("Unknown parameter "+arg);
				assert(false) : "Unknown parameter "+arg;
			}
		}
		assert(!inFiles.isEmpty()) : "No input vector files.";
		assert(!outNames.isEmpty()) : "No output (out=file[:frac][,file:frac...]).";
		assert(samplerate>0 && samplerate<=1) : "samplerate must be in (0,1]: "+samplerate;
		assert(balance>=0 && balance<1) : "balance must be in [0,1): "+balance;
		normalizeFractions();
	}

	/** Parse out=name[:frac],name[:frac],... — the last ':' introduces the fraction if it parses as a float. */
	private void parseOut(String b){
		for(String part : b.split(",")){
			if(part.length()==0){continue;}
			int colon=part.lastIndexOf(':');
			String name=part; float frac=-1;
			if(colon>0 && colon<part.length()-1){
				String maybe=part.substring(colon+1);
				Float f=tryFloat(maybe);
				if(f!=null){name=part.substring(0, colon); frac=f;}
			}
			outNames.add(name);
			outFracs.add(frac);
		}
	}

	private static Float tryFloat(String s){
		try{return Float.valueOf(s);}catch(NumberFormatException e){return null;}
	}

	/** Fill unspecified fractions and normalize the set to sum 1. */
	private void normalizeFractions(){
		int n=outFracs.size();
		float sumSpec=0; int nUnspec=0;
		for(float f : outFracs){if(f<0){nUnspec++;}else{sumSpec+=f;}}
		if(nUnspec>0){
			float rem=Math.max(0, 1f-sumSpec);
			float each=rem/nUnspec;
			for(int i=0; i<n; i++){if(outFracs.get(i)<0){outFracs.set(i, each);}}
		}
		float total=0;
		for(float f : outFracs){total+=f;}
		assert(total>0) : "Output fractions sum to 0.";
		for(int i=0; i<n; i++){outFracs.set(i, outFracs.get(i)/total);}
	}

	/*--------------------------------------------------------------*/

	void process(Timer t){
		Random randy=(seed>=0) ? new Random(seed) : new Random();

		//Read every input; capture the first #dims header; keep body rows only.
		ArrayList<byte[]> rows=new ArrayList<byte[]>();
		for(String f : inFiles){
			FileFormat ff=FileFormat.testInput(f, FileFormat.TXT, null, true, true);
			ByteFile bf=ByteFile.makeByteFile(ff);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='#'){
					if(dimsHeader==null && Tools.startsWith(line, "#dims")){dimsHeader=line.clone();}
					continue;//skip all header/comment lines
				}
				rows.add(line);
			}
			errorState|=bf.close();
		}
		long nIn=rows.size();

		if(dedupe){rows=deduplicate(rows);}
		if(samplerate<1f){rows=subsample(rows, samplerate, randy);}
		if(shuffle){Collections.shuffle(rows, randy);}

		//Partition into the output sets by (normalized) fraction.
		ArrayList<ArrayList<byte[]>> parts=partition(rows);

		long nOut=0;
		for(int i=0; i<outNames.size(); i++){
			ArrayList<byte[]> part=parts.get(i);
			if(balance>0){balance(part, randy);}
			if(balance>0 && shuffle){Collections.shuffle(part, randy);}//spread the duplicates
			write(outNames.get(i), part);
			int p=0; for(byte[] r : part){if(isPositive(r)){p++;}}
			outstream.println(outNames.get(i)+": "+part.size()+" rows ("+p+" pos, "+(part.size()-p)+" neg)");
			nOut+=part.size();
		}
		t.stop();
		outstream.println("Rows in: "+nIn+"  rows out: "+nOut+"  inputs: "+inFiles.size()+"  outputs: "+outNames.size());
		outstream.println("Time: \t"+t);
		if(errorState){throw new RuntimeException(getClass().getName()+" terminated in an error state.");}
	}

	/** Sort-based exact-duplicate removal. */
	private ArrayList<byte[]> deduplicate(ArrayList<byte[]> rows){
		byte[][] arr=rows.toArray(new byte[0][]);
		Arrays.sort(arr, BYTE_CMP);
		ArrayList<byte[]> out=new ArrayList<byte[]>(arr.length);
		for(int i=0; i<arr.length; i++){
			if(i==0 || BYTE_CMP.compare(arr[i], arr[i-1])!=0){out.add(arr[i]);}
		}
		return out;
	}

	private ArrayList<byte[]> subsample(ArrayList<byte[]> rows, float rate, Random randy){
		ArrayList<byte[]> out=new ArrayList<byte[]>((int)(rows.size()*rate)+16);
		for(byte[] r : rows){if(randy.nextFloat()<rate){out.add(r);}}
		return out;
	}

	private ArrayList<ArrayList<byte[]>> partition(ArrayList<byte[]> rows){
		int n=rows.size();
		ArrayList<ArrayList<byte[]>> parts=new ArrayList<ArrayList<byte[]>>();
		int idx=0; float cum=0;
		for(int i=0; i<outFracs.size(); i++){
			cum+=outFracs.get(i);
			int end=(i==outFracs.size()-1) ? n : Math.min(n, Math.round(cum*n));
			if(end<idx){end=idx;}
			parts.add(new ArrayList<byte[]>(rows.subList(idx, end)));
			idx=end;
		}
		return parts;
	}

	/** Upsample the minority class (label = last column >= 0.5) to fraction 'balance' of the total. */
	private void balance(ArrayList<byte[]> rows, Random randy){
		ArrayList<byte[]> pos=new ArrayList<byte[]>(), neg=new ArrayList<byte[]>();
		for(byte[] r : rows){(isPositive(r) ? pos : neg).add(r);}
		ArrayList<byte[]> minority=(pos.size()<=neg.size()) ? pos : neg;
		int min=minority.size(), maj=Math.max(pos.size(), neg.size());
		if(min==0 || maj==0){return;}//single-class: cannot balance
		int target=(int)(balance*maj/(1f-balance));//desired minority count
		for(int added=min; added<target; added++){
			rows.add(minority.get(randy.nextInt(min)));
		}
	}

	/** True if the row's last tab-delimited field is >= 0.5 (the positive label). */
	private boolean isPositive(byte[] row){
		lp.set(row);
		int terms=lp.terms();
		return terms>0 && lp.parseFloat(terms-1)>=0.5f;
	}

	private void write(String name, ArrayList<byte[]> rows){
		FileFormat ff=FileFormat.testOutput(name, FileFormat.TXT, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		ByteBuilder bb=new ByteBuilder(1<<16);
		if(dimsHeader!=null){bb.append(dimsHeader).nl();}
		for(byte[] r : rows){
			bb.append(r).nl();
			if(bb.length()>=(1<<16)){bsw.print(bb.toBytes()); bb.clear();}
		}
		if(bb.length()>0){bsw.print(bb.toBytes());}
		errorState|=bsw.poisonAndWait();
	}

	/*--------------------------------------------------------------*/

	private final ArrayList<String> inFiles=new ArrayList<String>();
	private final ArrayList<String> outNames=new ArrayList<String>();
	private final ArrayList<Float> outFracs=new ArrayList<Float>();
	private boolean shuffle=false, dedupe=false, overwrite=true;
	private float samplerate=1f, balance=0f;
	private long seed=1;
	private byte[] dimsHeader=null;

	private final LineParser1 lp=new LineParser1('\t');
	private PrintStream outstream=System.err;
	public boolean errorState=false;

	/** Lexicographic byte-array comparator for sort-based dedupe. */
	private static final java.util.Comparator<byte[]> BYTE_CMP=new java.util.Comparator<byte[]>(){
		public int compare(byte[] x, byte[] y){
			int n=Math.min(x.length, y.length);
			for(int i=0; i<n; i++){
				int d=(x[i]&0xFF)-(y[i]&0xFF);
				if(d!=0){return d;}
			}
			return x.length-y.length;
		}
	};

}
