package prot;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/** Canonical protein k-mer counts, per retained-member family or per query record.
 * Numeric encoding is native variable-width A48, not DDL fixed-triplet encoding.
 * Protein alphabet/strand geometry differs from DNA Clade spectra, so the schema
 * is explicitly tagged instead of pretending to be a Clade file.
 * @author Sayu
 */
public final class ProteinSpectrum{

	public static void main(String[] args){
		try{run(args);}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}
	public static void run(String[] args) throws Exception{
		final PreParser pp=new PreParser(args,ProteinSpectrum.class,false);
		try{new ProteinSpectrum(pp.args,pp.outstream).process();}finally{Shared.closeStream(pp.outstream);}
	}

	private ProteinSpectrum(String[] args,PrintStream stream) throws Exception{
		outstream=stream; parser.overwrite=false;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			final String key=pair[0].toLowerCase(Locale.ROOT),value=pair.length>1 ? pair[1] : null;
			if(key.equals("manifest")){manifest=value;}
			else if(key.equals("manifestsha80")){manifestPin=DigestSuffix.requireSuffix(value,key);}
			else if(key.equals("insha80")){inputPin=DigestSuffix.requireSuffix(value,key);}
			else if(key.equals("includeconsensus")){includeConsensus=Parse.parseBoolean(value);}
			else if(key.equals("namespace")){namespace=value;}
			else if(key.equals("k")){k=Parse.parseIntKMG(value);}
			else if(key.equals("shard")){shard=Parse.parseIntKMG(value);}
			else if(key.equals("shards")){shards=Parse.parseIntKMG(value);}
			else if(!parser.parse(arg,key,value)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		require((manifest!=null) != (parser.in1!=null),"Specify exactly one of manifest=MODELS_TSV or in=FASTA");
		require(parser.out1!=null && parser.out2==null && parser.in2==null,"Require one out=TSV; paired inputs/outputs unsupported");
		require(!parser.append,"Appending would duplicate schema/IDs; append unsupported");
		require(k>=1 && k<=3 && shards>0 && shard>=0 && shard<shards,"Require k=1..3 and 0<=shard<shards");
		require(manifest!=null || (shard==0 && shards==1),"Shards select family manifest rows, not query records");
		require(manifest==null || parser.maxReads<0,"Family spectra require all retained members; reads applies only to in=FASTA");
		require(manifestPin==null || manifest!=null,"manifestsha80 requires manifest");
		require(inputPin==null || parser.in1!=null,"insha80 requires in");
		require(!includeConsensus || (manifest!=null && manifestPin!=null),"includeconsensus requires manifest and manifestsha80");
		if(namespace==null){namespace=manifest==null ? "query" : "tail";}
		require(namespace.matches("[A-Za-z0-9_-]+"),"Namespace must be a safe token");
		require(Tools.testInputFiles(false,true,manifest,parser.in1),"Unreadable input");
		require(Tools.testOutputFiles(parser.overwrite,false,false,parser.out1),"Unwritable output");
		checkDistinct(manifest==null ? parser.in1 : manifest,parser.out1);
		Shared.AMINO_IN=true;
	}

	private void process() throws Exception{
		final Timer timer=new Timer();
		if(manifestPin!=null){verify(manifest,manifestPin);}
		if(inputPin!=null){verify(parser.in1,inputPin);}
		final ArrayList<Family> families=manifest==null ? null : families();
		final Spectrum[] spectra=families==null ? null : new Spectrum[families.size()];
		long memberRecords=0,memberResidues=0,consensusRecords=0,consensusResidues=0;
		if(families!=null){
			final AtomicInteger next=new AtomicInteger();
			final ArrayList<Thread> workers=new ArrayList<Thread>();
			for(int t=0;t<Math.min(parser.workers(),families.size());t++){
				workers.add(new Thread(new Runnable(){public void run(){
					try{
						for(int i=next.getAndIncrement();i<families.size();i=next.getAndIncrement()){
							final Family f=families.get(i);
							final Spectrum s=countFile(f.path,k);
							require(s.records==f.members && s.residues==f.residues,"Member/residue mismatch: "+f.id);
							verify(f.path,f.pin);
							if(includeConsensus){addConsensus(f,s);}
							spectra[i]=s;
						}
					}catch(Throwable failure){failure.printStackTrace(); KillSwitch.kill("ProteinSpectrum family worker failed");}
				}},"ProteinSpectrum-"+t));
			}
			ThreadWaiter.startAndWait(workers);
			if(manifestPin!=null){verify(manifest,manifestPin);}
			for(int i=0;i<families.size();i++){
				final Family f=families.get(i); final Spectrum s=spectra[i];
				require(s!=null && s.records==Math.addExact(f.members,includeConsensus ? 1 : 0),"Expected original members plus exactly one requested consensus: "+f.id);
				memberRecords=Math.addExact(memberRecords,f.members); memberResidues=Math.addExact(memberResidues,f.residues);
				consensusRecords=Math.addExact(consensusRecords,s.records-f.members);
				consensusResidues=Math.addExact(consensusResidues,s.residues-f.residues);
			}
		}
		final ByteStreamWriter writer=new ByteStreamWriter(parser.out1,parser.overwrite,false,true); writer.start();
		long rows=0,total=0,invalid=0;
		final ByteBuilder bb=new ByteBuilder();
		try{
			writer.println(header(k));
			if(families!=null){
				for(int i=0;i<families.size();i++){
					final Spectrum s=spectra[i]; require(s!=null,"Missing completed family spectrum");
					writeRow(bb,namespace,families.get(i).id,s); writer.print(bb);
					rows++; total=Math.addExact(total,s.total); invalid=Math.addExact(invalid,s.invalid);
				}
			}else{
				final Streamer reader=open(parser.in1,parser.maxReads,parser.threadsIn); reader.start();
				final HashSet<String> ids=new HashSet<String>(); final Spectrum s=new Spectrum(k);
				try{
					for(ListNum<Read> batch=reader.nextList();batch!=null;batch=reader.nextList()){
						for(Read r : batch){
							require(r!=null && r.bases!=null && r.mate==null,"Expected unpaired protein record");
							final String id=firstToken(r.id); require(ids.add(id),"Duplicate query ID: "+id);
							s.clear(); s.add(r.bases); writeRow(bb,namespace,id,s); writer.print(bb);
							rows++; total=Math.addExact(total,s.total); invalid=Math.addExact(invalid,s.invalid);
						}
					}
				}finally{require(!ReadWrite.closeStream(reader),"Protein input reader failed");}
				if(inputPin!=null){verify(parser.in1,inputPin);}
			}
		}finally{require(!writer.poisonAndWait(),"Spectrum output writer failed");}
		timer.stop(); outstream.println("Spectra: "+rows+"; valid k-mers: "+total+"; noncanonical residues: "+invalid+"; time: "+timer);
		if(families!=null){outstream.println("original_member_records="+memberRecords+" original_member_residues="+memberResidues+
			" added_consensus_records="+consensusRecords+" added_consensus_residues="+consensusResidues);}
	}

	/** Reads the actual HbmTailLibrary models.tsv contract; never infers membership from DDL registers. */
	private ArrayList<Family> families() throws Exception{
		final Path base=Paths.get(manifest).toAbsolutePath().getParent();
		final FamilyOutputGuard outputGuard=new FamilyOutputGuard(parser.out1);
		final ByteFile reader=ByteFile.makeByteFile(manifest,false);
		final ArrayList<Family> list=new ArrayList<Family>(); final HashSet<String> ids=new HashSet<String>();
		final LineParser1 lp=new LineParser1('\t');
		try{
			final byte[] first=reader.nextLine(); require(first!=null,"Empty family manifest");
			final String[] cols=new String(first,StandardCharsets.US_ASCII).split("\t",-1);
			final int idc=column(cols,"family_id"),pathc=column(cols,"directory"),mc=column(cols,"members"),rc=column(cols,"residues"),pc=column(cols,"retained_sha80");
			final int cc=includeConsensus ? column(cols,"consensus_sha80") : -1;
			int row=0;
			for(byte[] line=reader.nextLine();line!=null;line=reader.nextLine()){
				lp.set(line); require(lp.terms()==cols.length,"Family manifest field count mismatch");
				final String id=lp.parseString(idc); validId(id); require(ids.add(id),"Duplicate family: "+id);
				if(row++%shards!=shard){continue;}
				Path dir=Paths.get(lp.parseString(pathc)); if(!dir.isAbsolute()){dir=base.resolve(dir);}
				final String path=dir.resolve("retained.faa").toString(); outputGuard.check(path);
				final String consensusPath=includeConsensus ? dir.resolve("consensus.faa").toString() : null;
				if(includeConsensus){outputGuard.check(consensusPath);}
				final long members=lp.parseLong(mc),residues=lp.parseLong(rc);
				require(members>0 && residues>=0,"Invalid expected family totals: "+id);
				list.add(new Family(id,path,members,residues,DigestSuffix.requireSuffix(lp.parseString(pc),"retained_sha80"),
					consensusPath,includeConsensus ? DigestSuffix.requireSuffix(lp.parseString(cc),"consensus_sha80") : null));
			}
		}finally{require(!reader.close(),"Family manifest reader failed");}
		require(!list.isEmpty(),"No families selected"); return list;
	}

	public static Spectrum countFile(String path,int k){
		final Spectrum result=new Spectrum(k); final Streamer reader=open(path,-1,1); reader.start();
		try{
			for(ListNum<Read> batch=reader.nextList();batch!=null;batch=reader.nextList()){
				for(Read r : batch){require(r!=null && r.bases!=null && r.mate==null,"Expected unpaired family protein"); result.add(r.bases);}
			}
		}finally{require(!ReadWrite.closeStream(reader),"Family protein reader failed");}
		result.validate(); return result;
	}

	/** The library contract binds a singleton first-token family ID and file pin.
	 * Members were already counted/verified; add one copy even if its sequence is
	 * identical to a member. Each add() resets the rolling key at the record boundary.
	 */
	private static void addConsensus(Family family,Spectrum spectrum){
		assert(spectrum.records==family.members && spectrum.residues==family.residues) : "Verify original membership before adding the consensus";
		verify(family.consensusPath,family.consensusPin);
		final Streamer reader=open(family.consensusPath,-1,1); reader.start(); int records=0;
		try{
			for(ListNum<Read> batch=reader.nextList();batch!=null;batch=reader.nextList()){
				for(Read r : batch){
					require(records==0,"Consensus must contain exactly one record: "+family.id);
					require(r!=null && r.bases!=null && r.bases.length>0 && r.mate==null,"Consensus must be nonempty and unpaired: "+family.id);
					require(family.id.equals(firstToken(r.id)),"Consensus family ID mismatch: "+family.id);
					spectrum.add(r.bases); records++;
				}
			}
		}finally{require(!ReadWrite.closeStream(reader),"Consensus reader failed: "+family.id);}
		require(records==1,"Consensus must contain exactly one record: "+family.id);
		verify(family.consensusPath,family.consensusPin); spectrum.validate();
	}

	/** Reusable base20 counts: key(left...right), in AminoAcid.numberToAcid order. */
	public static final class Spectrum{
		public Spectrum(int order){require(order>=1 && order<=3,"Spectrum supports k=1..3"); k=order; counts=new long[POW[k]];}
		public void clear(){Arrays.fill(counts,0); total=records=residues=invalid=0;}
		public void add(byte[] seq){
			assert(seq!=null) : "Null is not a protein; each record resets the rolling k-mer";
			records=Math.addExact(records,1); residues=Math.addExact(residues,seq.length);
			int key=0,len=0;
			for(byte b : seq){
				final int aa=b>=0 && b<AminoAcid.acidToNumber.length ? AminoAcid.acidToNumber[b] : -1;
				if(aa<0 || aa>=20){invalid++; key=0; len=0; continue;}
				key=(key%POW[k-1])*20+aa; len=Math.min(k,len+1);
				if(len==k){counts[key]=Math.addExact(counts[key],1); total=Math.addExact(total,1);}
			}
		}
		public void validate(){
			long sum=0; for(long c : counts){require(c>=0,"Negative spectrum count"); sum=Math.addExact(sum,c);}
			require(total==sum,"Spectrum total must equal sum of all k-mer counts");
		}
		public final int k;
		public final long[] counts;
		public long total,records,residues,invalid;
	}

	/** Zero vectors return zero; callers retain a separate no-candidate sentinel. */
	public static double cosine(Spectrum a,Spectrum b){
		return cosine(a,b,norm(a),norm(b));
	}
	/** Precompute reference norms once; compute each query norm once before its shortlist. */
	public static double norm(Spectrum s){
		double squared=0; for(long count : s.counts){final double x=count; squared+=x*x;} return Math.sqrt(squared);
	}
	public static double cosine(Spectrum a,Spectrum b,double normA,double normB){
		require(a.k==b.k,"Cosine needs the same alphabet and k");
		require(Double.isFinite(normA) && Double.isFinite(normB) && normA>=0 && normB>=0,"Invalid cosine norms");
		if(normA==0 || normB==0){return 0;}
		double ab=0;
		for(int i=0;i<a.counts.length;i++){
			ab+=(double)a.counts[i]*b.counts[i];
		}
		return Math.min(1,ab/normA/normB);
	}

	public static String header(int k){return "#ProteinSpectrum\t1\tA48\tk="+k+"\talphabet="+ALPHABET;}
	public static void writeRow(ByteBuilder bb,String namespace,String id,Spectrum s){
		validId(id); require(namespace.matches("[A-Za-z0-9_-]+"),"Invalid namespace"); s.validate();
		bb.clear(); bb.append(namespace).tab().append(id).tab().appendA48(s.total).tab().appendA48(s.counts,'\t').nl();
	}

	/** Loads explicitly tagged protein spectra, preserving order. Keys are namespace:ID. */
	public static LinkedHashMap<String,Spectrum> load(String path){
		final ByteFile reader=ByteFile.makeByteFile(path,false); final LineParser1 lp=new LineParser1('\t');
		final LinkedHashMap<String,Spectrum> map=new LinkedHashMap<String,Spectrum>();
		try{
			final byte[] first=reader.nextLine(); require(first!=null,"Empty spectrum file");
			final String h=new String(first,StandardCharsets.US_ASCII);
			int order=0; for(int k=1;k<=3;k++){if(h.equals(header(k))){order=k;}}
			require(order>0,"Unsupported protein spectrum schema, alphabet, k or encoding");
			for(byte[] line=reader.nextLine();line!=null;line=reader.nextLine()){
				lp.set(line); require(lp.terms()==POW[order]+3,"Spectrum field count mismatch");
				final String ns=lp.parseString(0),id=lp.parseString(1); validId(id);
				require(ns.matches("[A-Za-z0-9_-]+"),"Invalid spectrum namespace");
				final Spectrum s=new Spectrum(order); s.total=a48(lp,2);
				for(int i=0;i<s.counts.length;i++){s.counts[i]=a48(lp,3+i);} s.validate();
				require(map.put(ns+":"+id,s)==null,"Duplicate spectrum: "+ns+":"+id);
			}
		}finally{require(!reader.close(),"Spectrum reader failed");}
		return map;
	}
	private static long a48(LineParser1 lp,int term){
		lp.setBounds(term); require(lp.b()>lp.a(),"Empty A48 count");
		long checked=0;
		for(int i=lp.a();i<lp.b();i++){
			final int digit=lp.line()[i]-48;
			require(digit>=0 && digit<64 && checked<=(Long.MAX_VALUE-digit)/64,"Invalid or overflowed nonnegative A48 count");
			checked=checked*64+digit;
		}
		final long nativeValue=lp.parseLongA48(term);
		assert(nativeValue==checked) : "Validated digits must roundtrip through native A48 decoding";
		return nativeValue;
	}
	private static Streamer open(String path,long limit,int threads){Shared.AMINO_IN=true; return StreamerFactory.makeStreamer(FileFormat.testInput(path,FileFormat.FASTA,null,true,false),null,true,limit,false,true,threads);}
	private static String firstToken(String id){
		require(id!=null,"Missing protein ID"); int end=0; while(end<id.length() && !Character.isWhitespace(id.charAt(end))){end++;}
		final String token=id.substring(0,end); validId(token); return token;
	}
	private static void validId(String id){require(id!=null && !id.isEmpty() && id.indexOf('\t')<0 && id.indexOf('\n')<0 && id.indexOf('\r')<0,"Missing/invalid spectrum ID");}
	private static int column(String[] names,String name){final int i=Arrays.asList(names).indexOf(name); require(i>=0,"Missing family manifest column: "+name); return i;}
	private static void verify(String path,String pin){require(DigestSuffix.file(path).equals(pin),"sha80 mismatch: "+path);}
	private static void checkDistinct(String a,String b) throws Exception{
		require(!new File(a).getCanonicalFile().equals(new File(b).getCanonicalFile()) &&
			!(Files.exists(Paths.get(a)) && Files.exists(Paths.get(b)) && Files.isSameFile(Paths.get(a),Paths.get(b))),"Input/output aliases: "+a+" and "+b);
	}
	/** Family mode finishes every input read, totals check and pin check before
	 * process() opens the output writer. A nonexistent output cannot alias any
	 * successfully read existing input; an unreadable input fails before writing.
	 * Keep the full alias checks for overwrite, including symlinks and hardlinks.
	 * This optimization does not apply to streaming per-record mode, whose writer
	 * opens before its reader. The global manifest/input guard stays unchanged.
	 */
	private static final class FamilyOutputGuard{
		FamilyOutputGuard(String name) throws Exception{
			assert(name!=null) : "Constructor requires the output checked by parameter validation";
			output=Paths.get(name); canonical=output.toFile().getCanonicalFile(); exists=Files.exists(output);
		}
		void check(String name) throws Exception{
			assert(name!=null) : "families() supplies each selected retained-member or consensus path";
			if(!exists){return;}
			final Path input=Paths.get(name);
			require(!input.toFile().getCanonicalFile().equals(canonical) &&
				!(Files.exists(input) && Files.isSameFile(input,output)),"Input/output aliases: "+name+" and "+output);
		}
		final Path output;
		final File canonical;
		final boolean exists;
	}
	private static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final class Family{
		Family(String id_,String path_,long members_,long residues_,String pin_,String consensusPath_,String consensusPin_){
			id=id_; path=path_; members=members_; residues=residues_; pin=pin_; consensusPath=consensusPath_; consensusPin=consensusPin_;
		}
		final String id,path,pin,consensusPath,consensusPin; final long members,residues;
	}
	public static final String ALPHABET=new String(AminoAcid.numberToAcid,0,20,StandardCharsets.US_ASCII);
	private static final int[] POW={1,20,400,8000};
	private final Parser parser=new Parser();
	private final PrintStream outstream;
	private String manifest,manifestPin,inputPin,namespace;
	private int k=2,shard=0,shards=1;
	private boolean includeConsensus=false;
}
