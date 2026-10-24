package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import map.IntLongHashMap;
import map.IntObjectMap;
import map.ObjectMap;
import map.ObjectSet;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/** Reports consensus-frame conservation and signed truth-end offsets.
 * No anchor, window width, or neural-network configuration is selected here.
 * Input records must already be sense-oriented, with actual lflank/rflank.
 * Unaligned reference overhangs are labeled extensions, not homologous columns.
 */
public final class NcrnaAnchorSlop {

	public static void main(String[] args){
		String fasta=null, consensus=null, models=null, prefix=null,exact=null,projectionCols=null;
		for(String arg:args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("fasta")||key.equals("in")){fasta=value;}
			else if(key.equals("consensus")){consensus=value;}
			else if(key.equals("models")){models=value;}
			else if(key.equals("outprefix")){prefix=value;}
			else if(key.equals("exact")){exact=value;}
			else if(key.equals("projectioncols")){projectionCols=value;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(fasta==null||consensus==null||models==null||prefix==null){
			throw new IllegalArgumentException("Required: fasta= consensus= models=name[,name] outprefix=; outputs must not exist");
		}
		run(fasta,consensus,models,prefix,exact,parseColumns(projectionCols));
	}

	public static void run(String fasta,String consensus,String requested,String prefix){
		run(fasta,consensus,requested,prefix,null);
	}

	static void run(String fasta,String consensus,String requested,String prefix,String exact){
		run(fasta,consensus,requested,prefix,exact,null);
	}

	public static int[] parseColumns(String text){
		if(text==null){return null;}
		final String[] fields=text.split(",",-1);
		final int[] columns=new int[fields.length];
		for(int i=0;i<fields.length;i++){
			columns[i]=Integer.parseInt(fields[i]);
			if(columns[i]<0){throw new IllegalArgumentException("Projection audit requires consensus columns >=0");}
			for(int j=0;j<i;j++){if(columns[j]==columns[i]){throw new IllegalArgumentException("Duplicate audit column: "+columns[i]);}}
		}
		return columns;
	}

	public static void run(String fasta,String consensus,String requested,String prefix,String exact,int[] auditColumns){
		for(String suffix:new String[]{".loci.tsv.gz",".columns.tsv",".offsets.tsv.gz",".counts.tsv"}){
			if(new File(prefix+suffix).exists()){throw new IllegalArgumentException("Output exists: "+prefix+suffix);}
		}
		if(auditColumns!=null){
			for(String suffix:new String[]{".projections.tsv.gz",".assignments.tsv.gz",".projection_summary.tsv"}){
				if(new File(prefix+suffix).exists()){throw new IllegalArgumentException("Output exists: "+prefix+suffix);}
			}
		}
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(consensus);
		final String[] loaded=TrnaConsensusBuilder.lastLoadedNames;
		if(library==null||library.length==0||loaded==null||loaded.length!=library.length){throw new IllegalArgumentException("Invalid consensus library");}
		final String[] names=new String[library.length];
		final HashSet<String> wanted=new HashSet<String>(Arrays.asList(requested.split(",",-1)));
		final Profile[] profiles=new Profile[library.length];
		for(int i=0;i<names.length;i++){
			names[i]=loaded[i].split("\\s+",2)[0];
			for(int j=0;j<i;j++){if(names[i].equals(names[j])){throw new IllegalArgumentException("Duplicate model: "+names[i]);}}
			if(wanted.remove(names[i])){profiles[i]=new Profile(names[i],library[i].length);}
		}
		if(!wanted.isEmpty()){throw new IllegalArgumentException("Unknown/empty requested models: "+wanted);}
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=true;
		final ArrayList<Read> reads=ReadInputStream.toReads(fasta,FileFormat.FA,-1);
		if(reads.isEmpty()){throw new IllegalArgumentException("Empty member FASTA");}
		final ObjectMap<String,byte[]> exactGenes=new ObjectMap<String,byte[]>(String.class,byte[].class);
		if(exact!=null){
			for(Read read:ReadInputStream.toReads(exact,FileFormat.FA,-1)){
				final String id=read.id.split("\\s+",2)[0];
				if(exactGenes.put(id,read.bases)!=null){throw new IllegalArgumentException("Duplicate exact gene ID: "+id);}
			}
			if(exactGenes.size()!=reads.size()){throw new IllegalArgumentException("Exact/flanked record counts differ");}
		}
		final ObjectSet<String> seenIds=new ObjectSet<String>(String.class);
		final ByteStreamWriter loci=writer(prefix+".loci.tsv.gz");
		final ProjectionAudit audit=auditColumns==null?null:new ProjectionAudit(prefix,auditColumns);
		loci.println("recordIndex\trecordId\tmodel\tassignment\tidentity\ttrueStart0\ttrueStop0\talignStart0\talignStop0\tstatus");
		long selected=0;
		try{
		for(int row=0;row<reads.size();row++){
			final Read read=reads.get(row);
			final int left=TrnaNinemerTableBuilder.parseFlankValue(read.id,"lflank="),right=TrnaNinemerTableBuilder.parseFlankValue(read.id,"rflank=");
			if(read.bases==null||left<0||right<0||(long)left+right>=read.length()){throw new IllegalArgumentException("Missing/invalid actual flanks: "+read.id);}
			final String recordId=read.id.split("\\s+",2)[0];
			if(!seenIds.add(recordId)){throw new IllegalArgumentException("Duplicate flanked gene ID: "+recordId);}
			if(exact!=null){requireExactGene(read.bases,left,right,exactGenes.get(recordId));}
			final String recorded=NcrnaConsensusPartitioner.recordedModel(read.id);
			int recordedIndex=-1;
			if(recorded!=null){
				for(int m=0;m<names.length;m++){if(names[m].equals(recorded)){recordedIndex=m;break;}}
				if(recordedIndex<0){throw new IllegalArgumentException("Unknown recorded model: "+read.id);}
			}
			final NcrnaConsensusPartitioner.Choice choice=NcrnaConsensusPartitioner.chooseModel(read.bases,library,recordedIndex);
			final Profile profile=profiles[choice.model];
			if(audit!=null){audit.assignment(row,recordId,names[choice.model],recorded!=null,profile!=null);}
			if(profile==null){continue;}
			selected++;profile.members++;
			final AlignmentStats stats=new AlignmentStats(true);
			final float identity=QuantumAligner.alignAndTraceStatic(library[choice.model],read.bases,stats);
			if(!Float.isFinite(identity)||stats.matchString==null){throw new IllegalStateException("No finite traced alignment: "+read.id);}
			final int[] projection=project(stats,library[choice.model].length,read.length());
			final int start=left,stop=read.length()-right-1;
			final boolean overlaps=stats.rStop>=start&&stats.rStart<=stop;
			if(audit!=null){audit.row(row,recordId,profile,read.bases,projection,start,stop,stats,identity,overlaps);}
			loci.println(row+"\t"+recordId+"\t"+profile.name+"\t"+(recorded==null?"highest_quantum_identity":"recorded_model")
				+"\t"+identity+"\t"+start+"\t"+stop+"\t"+stats.rStart+"\t"+stats.rStop+"\t"+(overlaps?"INCLUDED":"OFF_TRUTH_ALIGNMENT"));
			if(!overlaps){profile.offTruth++;continue;}
			for(int q=0;q<projection.length;q++){
				if(projection[q]>=0){profile.add(q,read.bases[projection[q]],start-projection[q],stop-projection[q]);}
			}
			// These positions are defined by the trace alone, without endpoint labels.
			for(int p=0;p<stats.rStart;p++){profile.add(p-stats.rStart,read.bases[p],start-p,stop-p);}
			for(int p=stats.rStop+1;p<read.length();p++){profile.add(profile.length+p-stats.rStop-1,read.bases[p],start-p,stop-p);}
		}
		}finally{if(audit==null){finish(loci);}else{finish(loci,audit.rows,audit.assignments);}}
		if(audit!=null){audit.summary(prefix,profiles);}
		final ByteStreamWriter columns=writer(prefix+".columns.tsv"),offsets=writer(prefix+".offsets.tsv.gz"),counts=writer(prefix+".counts.tsv");
		columns.println("model\tcolumn0\tframe\tmembers\tprojectable\tacgt\tprojectableFraction\tacgtCoverage\tmajorBase\tmajorFraction\tnormalizedEntropy\tstartMin\tstartMedian\tstartMax\tstopMin\tstopMedian\tstopMax");
		offsets.println("model\tcolumn0\tendpoint\tendpointMinusAnchor\tloci");
		counts.println("model\tmembers\toffTruthAlignment\tmodelLength");
		try{
		for(Profile profile:profiles){
			if(profile==null){continue;}
			counts.println(profile.name+"\t"+profile.members+"\t"+profile.offTruth+"\t"+profile.length);
			final int[] positions=profile.columns.toArray();Arrays.sort(positions);
			for(int position:positions){
				final Column c=profile.columns.get(position);
				final String frame=position<0?"LEFT_EXTENSION":position>=profile.length?"RIGHT_EXTENSION":"ALIGNED_COLUMN";
				columns.println(profile.name+"\t"+position+"\t"+frame+"\t"+profile.members+"\t"+c.projectable+"\t"+c.acgt()
					+"\t"+fraction(c.projectable,profile.members)+"\t"+fraction(c.acgt(),profile.members)+"\t"+c.majorBase()
					+"\t"+fraction(c.maximum(),c.acgt())+"\t"+c.entropy()+"\t"+range(c.start)+"\t"+range(c.stop));
				writeOffsets(offsets,profile.name,position,"START",c.start);
				writeOffsets(offsets,profile.name,position,"STOP",c.stop);
			}
		}
		}finally{finish(columns,offsets,counts);}
		System.err.println("ANCHOR_SLOP_DONE input="+reads.size()+" selected="+selected+" exactGeneParity="+(exact==null?"NOT_CHECKED":reads.size())+"; no anchor or NN width selected");
	}

	/** Additive diagnostic only: never changes assignment, alignment, or histogram inclusion. */
	public static final class ProjectionAudit {
		ProjectionAudit(String prefix,int[] columns_){
			if(columns_.length==0){throw new IllegalArgumentException("Empty audit columns");}
			columns=columns_;
			rows=writer(prefix+".projections.tsv.gz");assignments=writer(prefix+".assignments.tsv.gz");
			rows.println("recordIndex\trecordId\tmodel\tcolumn0\tstatus\tbase\tprojected0\ttrueStart0\ttrueStop0\tstartMinusAnchor\tstopMinusAnchor\talignStart0\talignStop0\toverlapBp\ttruthLength\talignedSpan\tidentity");
			assignments.println("recordIndex\trecordId\tmodel\tassignment\tprofileSelected");
		}
		void assignment(int row,String id,String model,boolean recorded,boolean selected){
			assert(row>=0&&id!=null&&model!=null):"Assignment audit must preserve the source record key";
			bb.clear();bb.append(row).tab().append(id).tab().append(model).tab().append(recorded?"recorded_model":"highest_quantum_identity").tab().append(selected?1:0).nl();assignments.print(bb);
		}
		void row(int row,String id,Profile profile,byte[] bases,int[] projection,int start,int stop,AlignmentStats stats,float identity,boolean overlaps){
			assert(projection.length==profile.length):"Audit projection must use the same model as the histogram";
			final int overlap=Math.max(0,Math.min(stop,stats.rStop)-Math.max(start,stats.rStart)+1);
			for(int col:columns){
				final int p=col<projection.length?projection[col]:-1;
				final byte base=p<0?(byte)'.':bases[p];
				final String status=!overlaps?"OFF_TRUTH_ALIGNMENT":col>=projection.length?"OUTSIDE_MODEL":p<0?"GAP":acgt(base)?"ACGT":"AMBIGUOUS";
				bb.clear();bb.append(row).tab().append(id).tab().append(profile.name).tab().append(col).tab().append(status).tab().append((char)base).tab();
				if(p<0){bb.append("NA");}else{bb.append(p);}
				bb.tab().append(start).tab().append(stop).tab();
				if(p<0){bb.append("NA\tNA");}else{bb.append(start-p).tab().append(stop-p);}
				bb.tab().append(stats.rStart).tab().append(stats.rStop).tab().append(overlap).tab().append(stop-start+1).tab().append(stats.rStop-stats.rStart+1).tab().append(identity,8).nl();rows.print(bb);
			}
		}
		void summary(String prefix,Profile[] profiles){
			final ByteStreamWriter out=writer(prefix+".projection_summary.tsv");
			out.println("model\tcolumn0\tmembers\toffTruthAlignment\tprojectable\tacgt\tunusable\tstartWithinMinus8Plus8\tstopWithinMinus8Plus8\tstartBelowMinus15\tstartAbovePlus9\tstopBelowMinus12");
			try{
				for(Profile profile:profiles){
					if(profile==null){continue;}
					for(int col:columns){
						final Column c=col<profile.length?profile.columns.get(col):null;
						final long n=c==null?0:c.acgt();
						bb.clear();bb.append(profile.name).tab().append(col).tab().append(profile.members).tab().append(profile.offTruth).tab().append(c==null?0:c.projectable).tab().append(n).tab().append(profile.members-n);
						bb.tab().append(c==null?0:mass(c.start,-8,8)).tab().append(c==null?0:mass(c.stop,-8,8));
						bb.tab().append(c==null?0:mass(c.start,Integer.MIN_VALUE,-16)).tab().append(c==null?0:mass(c.start,10,Integer.MAX_VALUE)).tab().append(c==null?0:mass(c.stop,Integer.MIN_VALUE,-13)).nl();out.print(bb);
					}
				}
			}finally{finish(out);}
		}
		public static long mass(IntLongHashMap histogram,int lo,int hi){
			assert(lo<=hi):"Inclusive histogram interval must not be reversed";
			long n=0;for(int key:histogram.toArray()){if(key>=lo&&key<=hi){n=Math.addExact(n,histogram.get(key));}}return n;
		}
		static boolean acgt(byte base){return base=='A'||base=='C'||base=='G'||base=='T';}
		final int[] columns;final ByteStreamWriter rows,assignments;final ByteBuilder bb=new ByteBuilder();
	}

	public static void requireExactGene(byte[] flanked,int left,int right,byte[] gene){
		if(gene==null||gene.length!=flanked.length-left-right){throw new IllegalArgumentException("Missing exact gene or de-flanked length mismatch");}
		for(int i=0;i<gene.length;i++){if(gene[i]!=flanked[left+i]){throw new IllegalArgumentException("De-flanked sequence differs from exact gene at offset "+i);}}
	}

	/** Consensus is query, sequence is reference. I consumes only query; D only reference. */
	public static int[] project(AlignmentStats stats,int queryLength,int sequenceLength){
		if(stats.matchString==null||stats.rStart<0||stats.rStop<stats.rStart||stats.rStop>=sequenceLength){throw new IllegalArgumentException("Invalid trace bounds");}
		final int[] map=new int[queryLength];Arrays.fill(map,-1);
		int q=0,r=stats.rStart;
		for(byte op:stats.matchString){
			if(op=='D'){r++;}
			else if(op=='I'){q++;}
			else if(op=='m'||op=='S'||op=='N'){
				if(q>=queryLength||r>stats.rStop){throw new IllegalArgumentException("Trace exceeds declared bounds");}
				map[q++]=r++;
			}else{throw new IllegalArgumentException("Unknown trace operation: "+(char)op);}
			if(q>queryLength||r>stats.rStop+1){throw new IllegalArgumentException("Trace consumption overflow");}
		}
		if(q!=queryLength||r!=stats.rStop+1){throw new IllegalArgumentException("Trace consumption mismatch q="+q+" r="+r);}
		return map;
	}

	static final class Profile {
		Profile(String name_,int length_){name=name_;length=length_;for(int i=0;i<length;i++){columns.put(i,new Column());}}
		void add(int column,byte base,int startOffset,int stopOffset){Column c=columns.get(column);if(c==null){c=new Column();columns.put(column,c);}c.add(base,startOffset,stopOffset);}
		final String name;final int length;long members,offTruth;
		final IntObjectMap<Column> columns=new IntObjectMap<Column>(Column.class);
	}
	public static final class Column {
		public void add(byte base,int startOffset,int stopOffset){projectable++;final int x=base=='A'?0:base=='C'?1:base=='G'?2:base=='T'?3:-1;if(x<0){return;}bases[x]++;start.increment(startOffset,1);stop.increment(stopOffset,1);}
		public long acgt(){return bases[0]+bases[1]+bases[2]+bases[3];}
		long maximum(){long max=0;for(long n:bases){max=Math.max(max,n);}return max;}
		public char majorBase(){if(acgt()==0){return '.';}for(int i=0;i<4;i++){if(bases[i]==maximum()){return "ACGT".charAt(i);}}throw new AssertionError();}
		public String entropy(){if(acgt()==0){return "NA";}double entropy=0;for(long n:bases){if(n>0){double p=(double)n/acgt();entropy-=p*Math.log(p)/Math.log(4);}}return String.format(Locale.ROOT,"%.8f",entropy);}
		public long projectable;final long[] bases=new long[4];
		public final IntLongHashMap start=new IntLongHashMap(),stop=new IntLongHashMap();
	}
	public static String range(IntLongHashMap histogram){
		if(histogram.isEmpty()){return "NA\tNA\tNA";}
		final int[] keys=histogram.toArray();Arrays.sort(keys);
		final long n=histogram.sum();long cumulative=0;int median=keys[0];
		for(int key:keys){cumulative+=histogram.get(key);if(cumulative>(n-1)/2){median=key;break;}}
		return keys[0]+"\t"+median+"\t"+keys[keys.length-1];
	}
	private static String fraction(long n,long total){return total==0?"NA":String.format(Locale.ROOT,"%.8f",(double)n/total);}
	private static void writeOffsets(ByteStreamWriter writer,String model,int column,String end,IntLongHashMap histogram){final int[] keys=histogram.toArray();Arrays.sort(keys);for(int key:keys){writer.println(model+"\t"+column+"\t"+end+"\t"+key+"\t"+histogram.get(key));}}
	private static ByteStreamWriter writer(String path){if(new File(path).exists()){throw new IllegalArgumentException("Output exists: "+path);}ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));w.start();return w;}
	private static void finish(ByteStreamWriter... writers){boolean error=false;for(ByteStreamWriter writer:writers){error|=writer.poisonAndWait();}if(error){throw new IllegalStateException("Diagnostic output write failed");}}
}
