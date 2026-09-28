package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import idaligner.QuantumAligner;
import map.ObjectIntMap;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/** Assigns flanked ncRNA genes to consensus models and writes one FASTA per model. */
public final class NcrnaConsensusPartitioner {

	public static void main(String[] args){
		String fasta=null,consensus=null,outDir=null,assignments=null,counts=null;
		for(String arg:args){final int eq=arg.indexOf('=');if(eq<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),value=arg.substring(eq+1);if(key.equals("fasta")||key.equals("in")){fasta=value;}else if(key.equals("consensus")){consensus=value;}else if(key.equals("outdir")){outDir=value;}else if(key.equals("assignments")||key.equals("assignmentsout")){assignments=value;}else if(key.equals("counts")||key.equals("countsout")){counts=value;}else{throw new IllegalArgumentException("Unknown argument: "+arg);}}
		if(fasta==null||consensus==null||outDir==null||assignments==null||counts==null){throw new IllegalArgumentException("Required: fasta= consensus= outdir= assignments= counts=");}
		new NcrnaConsensusPartitioner(fasta,consensus,outDir,assignments,counts).run();
	}

	NcrnaConsensusPartitioner(String fasta_,String consensus_,String outDir_,String assignments_,String counts_){fasta=fasta_;consensus=consensus_;outDir=outDir_;assignments=assignments_;counts=counts_;}

	void run(){
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(consensus);final String[] loaded=TrnaConsensusBuilder.lastLoadedNames;
		if(library==null||library.length<1||loaded==null||loaded.length!=library.length){throw new IllegalArgumentException("Consensus library/name mismatch: "+consensus);}
		final String[] names=new String[loaded.length],files=new String[loaded.length];final ObjectIntMap<String> nameToIndex=new ObjectIntMap<String>(String.class);
		for(int i=0;i<names.length;i++){names[i]=firstToken(loaded[i]);files[i]=safeName(names[i]);if(nameToIndex.put(names[i],i)>=0){throw new IllegalArgumentException("Duplicate consensus name: "+names[i]);}for(int j=0;j<i;j++){if(files[i].equals(files[j])){throw new IllegalArgumentException("Consensus names collide as filenames: "+names[j]+" and "+names[i]);}}}
		final File dir=new File(outDir);if(!dir.isDirectory()&&!dir.mkdirs()){throw new IllegalArgumentException("Could not create output directory: "+outDir);}
		final ByteStreamWriter[] writers=new ByteStreamWriter[names.length];for(int i=0;i<writers.length;i++){writers[i]=writer(new File(dir,files[i]+".fa.gz").getPath());writers[i].start();}
		final ByteStreamWriter assignmentWriter=writer(assignments);assignmentWriter.start();final ByteBuilder assignmentBuffer=new ByteBuilder(1<<16).append("recordIndex\trecordId\tmodelIndex\tmodel\tassignmentSource\tidentity\n");
		Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;final ArrayList<Read> reads=ReadInputStream.toReads(fasta,FileFormat.FA,-1);if(reads.isEmpty()){throw new IllegalArgumentException("Empty flanked FASTA: "+fasta);}
		final long[] modelCounts=new long[names.length];long recorded=0,fallback=0;
		for(int record=0;record<reads.size();record++){
			final Read read=reads.get(record);validateFlanks(read);
			final String recordedName=recordedModel(read.id);final int recordedIndex=(recordedName==null ? -1 : nameToIndex.get(recordedName));
			if(recordedName!=null&&recordedIndex<0){throw new IllegalArgumentException("Recorded accepting model is absent from consensus library: "+recordedName+" header="+read.id);}
			final Choice choice=chooseModel(read.bases,library,recordedIndex);final String source=(recordedIndex>=0 ? "recorded_accepting_model" : "highest_quantum_identity");if(recordedIndex>=0){recorded++;}else{fallback++;}
			modelCounts[choice.model]++;final ByteBuilder fastaBuffer=new ByteBuilder(read.length()+read.id.length()+128).append('>').append(read.id);if(tokenValue(read.id,"accepted_model=")==null){fastaBuffer.append(" accepted_model=").append(names[choice.model]);}fastaBuffer.append(" model_assignment_source=").append(source).append(" model_identity=").append(choice.identity,9).nl().append(read.bases).nl();writers[choice.model].print(fastaBuffer);
			assignmentBuffer.append(record).tab().append(firstToken(read.id)).tab().append(choice.model).tab().append(names[choice.model]).tab().append(source).tab().append(choice.identity,9).nl();if(assignmentBuffer.length>=1<<15){assignmentWriter.print(assignmentBuffer);assignmentBuffer.clear();}
		}
		assignmentWriter.print(assignmentBuffer);if(assignmentWriter.poisonAndWait()){throw new RuntimeException("Write error: "+assignments);}for(ByteStreamWriter out:writers){if(out.poisonAndWait()){throw new RuntimeException("Partition FASTA write error");}}
		final ByteStreamWriter countWriter=writer(counts);countWriter.start();final ByteBuilder countBuffer=new ByteBuilder().append("modelIndex\tmodel\trecords\ttrainingEligible500\n");long total=0;for(int i=0;i<names.length;i++){countBuffer.append(i).tab().append(names[i]).tab().append(modelCounts[i]).tab().append(modelCounts[i]>=500).nl();total+=modelCounts[i];}countBuffer.append("#total\t.\t").append(total).tab().append('.').nl().append("#recorded\t.\t").append(recorded).tab().append('.').nl().append("#highest_identity_fallback\t.\t").append(fallback).tab().append('.').nl();countWriter.print(countBuffer);if(countWriter.poisonAndWait()){throw new RuntimeException("Write error: "+counts);}
		if(total!=reads.size()){throw new IllegalStateException("Assignment conservation failed: total="+total+" reads="+reads.size());}
		System.err.println("NCRNA_CONSENSUS_PARTITION_PASS records="+total+" recorded="+recorded+" highestIdentityFallback="+fallback);
	}

	static Choice chooseModel(byte[] bases,byte[][] library,int recordedModel){
		if(recordedModel>=library.length){throw new IllegalArgumentException("Recorded model index out of range: "+recordedModel);}
		if(recordedModel>=0){return new Choice(recordedModel,QuantumAligner.alignStatic(library[recordedModel],bases,null));}
		int best=-1;float bestIdentity=-1;for(int i=0;i<library.length;i++){final float identity=QuantumAligner.alignStatic(library[i],bases,null);if(identity>bestIdentity){best=i;bestIdentity=identity;}}
		if(best<0||!Float.isFinite(bestIdentity)){throw new IllegalStateException("No finite consensus identity");}return new Choice(best,bestIdentity);
	}

	static String recordedModel(String header){String value=tokenValue(header,"accepted_model=");if(value==null){value=tokenValue(header,"model=");}return value;}
	private static String tokenValue(String header,String key){int from=0;while(from<header.length()){final int idx=header.indexOf(key,from);if(idx<0){return null;}if(idx==0||Character.isWhitespace(header.charAt(idx-1))){final int start=idx+key.length();int stop=start;while(stop<header.length()&&!Character.isWhitespace(header.charAt(stop))){stop++;}if(stop>start){return header.substring(start,stop);}}from=idx+key.length();}return null;}
	private static void validateFlanks(Read read){final int left=TrnaNinemerTableBuilder.parseFlankValue(read.id,"lflank="),right=TrnaNinemerTableBuilder.parseFlankValue(read.id,"rflank=");if(left<0||right<0||left+right>=read.length()){throw new IllegalArgumentException("Invalid flanks in header: "+read.id);}}
	private static String safeName(String name){final StringBuilder sb=new StringBuilder(name.length());for(int i=0;i<name.length();i++){final char c=name.charAt(i);sb.append(Character.isLetterOrDigit(c)||c=='.'||c=='_'||c=='-' ? c : '_');}if(sb.length()<1){throw new IllegalArgumentException("Empty consensus filename for "+name);}return sb.toString();}
	private static String firstToken(String s){int i=0;while(i<s.length()&&!Character.isWhitespace(s.charAt(i))){i++;}return s.substring(0,i);}
	private static ByteStreamWriter writer(String path){return new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));}

	static final class Choice{Choice(int model_,float identity_){model=model_;identity=identity_;}final int model;final float identity;}
	final String fasta,consensus,outDir,assignments,counts;
}
