package prot;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import map.LongIntMap;
import map.LongLongHashMap;
import map.ObjectIntMap;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Builds the per-reference-CDS contribution table used by simulated breakpoint
 * clones. The GFF stream and per-tid gene ordinal assignment come directly from
 * {@link ReferenceCdsShredSurvivalTable#scanGenes}; this class adds the exact
 * CallGenes protein identity and canonical family assignment for each physical
 * CDS row.
 *
 * <p>Each output row is keyed by {@code (tid,gene_ordinal)} and retains the
 * physical CDS contribution that must be removed with that gene: one CDS,
 * {@code length}, {@code length^2}, coding length, and an optional assigned
 * family rank. Multipart identities therefore occupy several rows with one
 * ordinal rather than being reduced to an interval-envelope approximation.</p>
 *
 * <p>{@code proteins=} is consumed in lockstep with the GFF. Every expected
 * {@code contig_gN} header must match exactly, so a contig-local index drift,
 * permutation, missing protein, or extra protein fails before publication.
	 * {@code assignmentmanifest=} binds both the sparse gene-hit file and complete
	 * assignment-provenance file for every chunk by sha80. Every query must appear
	 * exactly once in the provenance, including masked/bait interceptions and typed
	 * rejections; only {@code ASSIGNED} rows carry a decrementable family rank.</p>
 *
 * @author Yoimiya
 */
public final class ReferenceCdsGeneContributionTable {

	static final String SCHEMA="reference_cds_gene_contribution_v1";
	static final String ASSIGNMENT_MANIFEST_SCHEMA="reference_cds_gene_assignment_manifest_v1";
	static final String COLUMNS="G\tquery_id\ttid\tgene_ordinal\tcontig_id\tlocal_gene_index\tstart0\tend0\tlength\tfamily_rank\tassignment_outcome";

	private ReferenceCdsGeneContributionTable(){}

	/** Summary returned by {@link #build}. */
	public static final class Summary {
		public long rows, identities, mapped, masked, bait, proteinHeaders, assignmentRows, hitRows;
		public int tids, contigs;
	}

	/** Command-line entry point. */
	public static void main(String[] args){
		String gff=null, assembly=null, proteins=null, assignmentManifest=null, roster=null, out=null;
		int tidFallback=-1, expectedFamilies=-1; boolean overwrite=false;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new RuntimeException("Arguments must be key=value: "+arg);}
			final String key=arg.substring(0,eq).toLowerCase(), value=arg.substring(eq+1);
			if(key.equals("gff")){gff=value;}
			else if(key.equals("assembly")){assembly=value;}
			else if(key.equals("proteins")){proteins=value;}
			else if(key.equals("assignmentmanifest")){assignmentManifest=value;}
			else if(key.equals("roster")){roster=value;}
			else if(key.equals("out")){out=value;}
			else if(key.equals("tid")){tidFallback=Integer.parseInt(value);}
			else if(key.equals("expectedfamilies")){expectedFamilies=Integer.parseInt(value);}
			else if(key.equals("overwrite") || key.equals("ow")){overwrite=parse.Parse.parseBoolean(value);}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(gff==null || assembly==null || proteins==null || assignmentManifest==null || roster==null || out==null || expectedFamilies<1){
			throw new RuntimeException("Usage: java -ea prot.ReferenceCdsGeneContributionTable gff=<whole.gff[.gz]> assembly=<whole.fa[.gz]> proteins=<geneid.faa[.gz]> assignmentmanifest=<assignments.tsv> roster=<roster.tsv> expectedfamilies=<n> out=<table.tsv> [tid=<fallback>] [ow=f]");
		}
		final Summary s=build(gff,assembly,proteins,assignmentManifest,roster,out,tidFallback,expectedFamilies,overwrite);
		System.err.println("ReferenceCdsGeneContributionTable PASS: rows="+s.rows+" identities="+s.identities
			+" mapped="+s.mapped+" masked="+s.masked+" bait="+s.bait+" tids="+s.tids+" contigs="+s.contigs
			+" proteins="+s.proteinHeaders+" assignments="+s.assignmentRows+" hits="+s.hitRows+" out="+out);
	}

	/** Builds one failure-safe contribution table. */
	public static Summary build(final String gffPath, final String assemblyPath, final String proteinsPath,
			final String assignmentManifestPath, final String rosterPath, final String outPath,
			final int tidFallback, final int expectedFamilies, final boolean overwrite){
		final File outFile=new File(outPath);
		for(String in : new String[]{gffPath,assemblyPath,proteinsPath,assignmentManifestPath,rosterPath}){
			if(ReferenceCdsShredSurvivalTable.samePath(in,outPath)){
				throw new IllegalArgumentException("out= names an input file: "+outPath);
			}
		}
		if(outFile.exists() && !overwrite){throw new IllegalArgumentException("out= already exists (pass ow=t to replace it): "+outPath);}
		final File partial;
		try{
			final File dir=outFile.getAbsoluteFile().getParentFile();
			if(dir!=null && !dir.isDirectory()){throw new IllegalArgumentException("Output directory does not exist: "+dir);}
			partial=File.createTempFile(outFile.getName()+".",".partial",dir);
		}catch(IOException e){throw new RuntimeException("Could not create temporary output beside "+outPath,e);}
		boolean done=false;
		try{
			final Summary summary=buildTo(gffPath,assemblyPath,proteinsPath,assignmentManifestPath,
				rosterPath,partial.getPath(),tidFallback,expectedFamilies);
			ReferenceCdsShredSurvivalTable.publishCompleted(partial,outFile,overwrite);
			done=true;
			return summary;
		}finally{
			if(!done && partial.exists()){partial.delete();}
		}
	}

	private static Summary buildTo(final String gffPath, final String assemblyPath, final String proteinsPath,
			final String assignmentManifestPath, final String rosterPath, final String outPath,
			final int tidFallback, final int expectedFamilies){
		final HashRankMap roster=loadRoster(rosterPath,expectedFamilies);
		final AssignmentIndex assignments=loadAssignments(assignmentManifestPath,roster);
		final ProteinContigIndex proteins=ProteinContigIndex.load(proteinsPath);
		final ByteStreamWriter writer=new ByteStreamWriter(outPath,true,false,true); writer.start();
		final ByteBuilder header=new ByteBuilder(2048);
		header.append("#schema_version\t").append(SCHEMA).nl();
		header.append("#tool\tprot.ReferenceCdsGeneContributionTable").nl();
		header.append("#ordinal_source\tprot.ReferenceCdsShredSurvivalTable.scanGenes").nl();
		header.append("#contribution_rule\tone row per physical GFF CDS; multipart identities share gene_ordinal; removing an ordinal removes every row carrying it").nl();
		header.append("#gff_sha80\t").append(gffPath).tab().append(sha80(gffPath)).nl();
		header.append("#assembly_sha80\t").append(assemblyPath).tab().append(sha80(assemblyPath)).nl();
		header.append("#proteins_sha80\t").append(proteinsPath).tab().append(sha80(proteinsPath)).nl();
		header.append("#assignment_manifest_sha80\t").append(assignmentManifestPath).tab().append(sha80(assignmentManifestPath)).nl();
		header.append("#roster_sha80\t").append(rosterPath).tab().append(sha80(rosterPath)).nl();
		header.append("#expected_families\t").append(expectedFamilies).nl();
		header.append("#columns\t").append(COLUMNS).nl();
		writer.print(header);
		final ContributionSink sink=new ContributionSink(writer,proteins,assignments);
		ReferenceCdsShredSurvivalTable.Summary genes=null;
		boolean ok=false;
		try{
			genes=ReferenceCdsShredSurvivalTable.scanGenes(gffPath,assemblyPath,tidFallback,sink);
			sink.finish();
			if(sink.rows!=genes.cdsRows){throw new AssertionError("Contribution rows "+sink.rows+" != canonical CDS rows "+genes.cdsRows);}
			if(assignments.consumed.size()!=assignments.ranks.size()){
				throw new IllegalArgumentException("Only "+assignments.consumed.size()+" of "+assignments.ranks.size()+" assignment-provenance queries matched the GFF/protein stream");
			}
			sink.flush();
			final ByteBuilder trailer=new ByteBuilder(256);
			trailer.append("#end\tG=").append(sink.rows).append("\tidentities=").append(genes.genes)
				.append("\tmapped=").append(sink.mapped).append("\tmasked=").append(sink.masked)
				.append("\tbait=").append(sink.bait).append("\ttids=").append(genes.tids)
				.append("\tcontigs=").append(genes.contigs).append("\tproteins=").append(proteins.headers)
				.append("\tassignments=").append(assignments.ranks.size()).append("\thits=").append(assignments.hitRows).nl();
			writer.print(trailer);
			ok=true;
		}finally{if(!ok){writer.poisonAndWait();}}
		if(writer.poisonAndWait()){throw new RuntimeException("I/O error writing "+outPath);}
		final Summary out=new Summary();
		out.rows=sink.rows; out.identities=genes.genes; out.mapped=sink.mapped; out.masked=sink.masked; out.bait=sink.bait;
		out.proteinHeaders=proteins.headers; out.assignmentRows=assignments.ranks.size(); out.hitRows=assignments.hitRows;
		out.tids=genes.tids; out.contigs=genes.contigs;
		return out;
	}

	/** Emits one physical CDS contribution without allocating a query-id String. */
	private static final class ContributionSink implements ReferenceCdsShredSurvivalTable.GeneContributionSink {
		final ByteStreamWriter writer; final ProteinContigIndex proteins; final AssignmentIndex assignments;
		final ByteBuilder query=new ByteBuilder(256), buffer=new ByteBuilder(1<<16);
		String currentContig=null; int currentGenes=0;
		long rows=0, mapped=0, masked=0, bait=0;
		ContributionSink(ByteStreamWriter writer_, ProteinContigIndex proteins_, AssignmentIndex assignments_){writer=writer_; proteins=proteins_; assignments=assignments_;}
		@Override
		public void observe(final String contig, final int localGeneIndex, final int tid,
				final int geneOrdinal, final int start0, final int end0){
			if(currentContig==null || !currentContig.equals(contig)){
				finishContig(); currentContig=contig; currentGenes=0;
			}
			if(localGeneIndex!=currentGenes){throw new IllegalArgumentException("Noncontiguous GFF local gene index for "+contig+": expected "+currentGenes+" got "+localGeneIndex);}
			currentGenes++;
			query.clear(); query.append(contig).append("_g").append(localGeneIndex);
			final int encoded=assignments.get(query.array,0,query.length);
			final Outcome outcome=Outcome.fromEncoded(encoded);
			final int family=outcome==Outcome.ASSIGNED ? (encoded&0xFFFF)-1 : -1;
			if(outcome==Outcome.ASSIGNED){mapped++;}
			else if(outcome==Outcome.MASKED_INTERCEPTED){masked++;}
			else if(outcome==Outcome.BAIT_INTERCEPTED){bait++;}
			final int length=end0-start0;
			assert(length>0) : "Nonpositive CDS length for "+contig+"_g"+localGeneIndex+": "+start0+"-"+end0;
			buffer.append('G').tab().append(query.array,0,query.length).tab().append(tid).tab().append(geneOrdinal).tab()
				.append(contig).tab().append(localGeneIndex).tab().append(start0).tab().append(end0).tab().append(length).tab();
			if(family<0){buffer.append('-');}else{buffer.append(family);}
			buffer.tab().append(outcome.text).nl(); rows++;
			if(buffer.length>=(1<<16)){writer.print(buffer); buffer.clear();}
		}
		void finish(){finishContig(); if(proteins.counts.size()!=0){throw new IllegalArgumentException(proteins.counts.size()+" protein contig(s) were absent from the GFF");} flush();}
		void flush(){if(buffer.length>0){writer.print(buffer); buffer.clear();}}
		private void finishContig(){
			if(currentContig==null){return;}
			final int expected=proteins.counts.get(currentContig);
			if(expected<0){throw new IllegalArgumentException("GFF contig absent from protein FASTA: "+currentContig);}
			if(expected!=currentGenes){throw new IllegalArgumentException("GFF/protein gene-count mismatch for "+currentContig+": GFF="+currentGenes+" protein="+expected);}
			proteins.counts.remove(currentContig); currentContig=null; currentGenes=0;
		}
	}

	/** Validates per-contig {@code _gN} header runs without assuming global contig order matches the GFF. */
	private static final class ProteinContigIndex {
		final ObjectIntMap<String> counts=new ObjectIntMap<String>(1<<16,String.class); long headers=0;
		static ProteinContigIndex load(final String path){
			final ProteinContigIndex out=new ProteinContigIndex();
			final ByteFile input=ByteFile.makeByteFile(path,true);
			byte[] current=null; int currentLength=0,currentCount=0; long lineNo=0; boolean ok=false;
			try{
				for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
					lineNo++; if(line.length==0 || line[0]!='>'){continue;}
					int end=line.length; while(end>1 && line[end-1]<=' '){end--;}
					int marker=end-2; while(marker>1 && line[marker]!='_'){marker--;}
					if(marker<=1 || marker+2>=end || line[marker+1]!='g'){throw bad(path,lineNo,"protein header lacks terminal _gN: "+new String(line));}
					final int index=parseHeaderIndex(line,marker+2,end,path,lineNo);
					final int contigLength=marker-1;
					boolean same=current!=null && currentLength==contigLength;
					if(same){for(int i=0; i<contigLength; i++){if(current[i]!=line[i+1]){same=false;break;}}}
					if(!same){
						if(current!=null){putProteinContig(out,current,currentLength,currentCount,path);}
						current=new byte[contigLength]; System.arraycopy(line,1,current,0,contigLength); currentLength=contigLength; currentCount=0;
					}
					if(index!=currentCount){throw bad(path,lineNo,"protein local index "+index+" != expected "+currentCount);}
					currentCount++; out.headers++;
				}
				if(current!=null){putProteinContig(out,current,currentLength,currentCount,path);}
				ok=true;
			}finally{if(input.close() && ok){throw new RuntimeException("I/O error reading "+path);}}
			if(out.headers==0){throw new IllegalArgumentException("Protein FASTA has no headers: "+path);}
			return out;
		}
		private static void putProteinContig(ProteinContigIndex out,byte[] key,int length,int count,String path){
			final String contig=new String(key,0,length);
			if(out.counts.put(contig,count)>=0){throw new IllegalArgumentException("Noncontiguous or duplicate protein contig in "+path+": "+contig);}
		}
		private static int parseHeaderIndex(byte[] line,int from,int to,String path,long lineNo){
			long value=0; if(from>=to){throw bad(path,lineNo,"empty _g index");}
			for(int i=from; i<to; i++){final byte c=line[i]; if(c<'0' || c>'9'){throw bad(path,lineNo,"non-digit in _g index");} value=value*10+(c-'0'); if(value>Integer.MAX_VALUE){throw bad(path,lineNo,"_g index overflow");}}
			return (int)value;
		}
	}

	/** Two independent 64-bit hashes make byte-key collisions detectable. */
	private static class HashRankMap {
		final LongIntMap ranks; final LongLongHashMap second;
		HashRankMap(int expected){ranks=new LongIntMap(Math.max(16,expected)); second=new LongLongHashMap(Math.max(16,expected));}
		void put(final byte[] key, final int from, final int to, final int rank, final String source){
			final long h1=hash1(key,from,to), h2=hash2(key,from,to);
			if(ranks.contains(h1)){
				if(second.get(h1)!=h2){throw new IllegalArgumentException("Primary query hash collision in "+source);}
				throw new IllegalArgumentException("Duplicate key in "+source+": "+new String(key,from,to-from));
			}
			ranks.put(h1,rank+1);
			if(!second.put(h1,h2)){throw new AssertionError("Secondary hash map rejected a new primary key");}
		}
		int get(final byte[] key, final int from, final int to){
			final long h1=hash1(key,from,to);
			final int value=ranks.get(h1);
			if(value<0){return -1;}
			final long h2=hash2(key,from,to);
			if(second.get(h1)!=h2){throw new IllegalArgumentException("Primary query hash collision during lookup");}
			return value-1;
		}
		int size(){return ranks.size();}
	}

	private static final class AssignmentIndex extends HashRankMap {
		final LongHashSet consumed;
		long hitRows=0;
		AssignmentIndex(int expected){super(expected); consumed=new LongHashSet(Math.max(16,expected));}
		int peek(final byte[] key, final int from, final int to){return super.get(key,from,to);}
		@Override
		int get(final byte[] key, final int from, final int to){
			final int encoded=super.get(key,from,to);
			if(encoded<0){throw new IllegalArgumentException("Query absent from assignment provenance: "+new String(key,from,to-from));}
			final long h1=hash1(key,from,to);
			if(!consumed.add(h1)){throw new IllegalArgumentException("Assignment query appeared more than once in the GFF/protein stream: "+new String(key,from,to-from));}
			return encoded;
		}
	}

	private static HashRankMap loadRoster(final String path, final int expectedFamilies){
		final HashRankMap out=new HashRankMap(expectedFamilies);
		final ByteFile input=ByteFile.makeByteFile(path,true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		long lineNo=0; boolean header=false, ok=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lineNo++;
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(!header){
					if(lp.terms()<3 || !lp.termEquals("active_index",0) || !lp.termEquals("rep_id",2)){
						throw bad(path,lineNo,"expected active_index/family_id/rep_id roster header");
					}
					header=true; continue;
				}
				if(lp.terms()<3){throw bad(path,lineNo,"roster row has fewer than three fields");}
				final int rank=lp.parseInt(0);
				if(rank!=out.size()){throw bad(path,lineNo,"active_index "+rank+" != contiguous row "+out.size());}
				lp.setBounds(2); out.put(line,lp.a(),lp.b(),rank,path+":"+lineNo);
			}
			ok=true;
		}finally{if(input.close() && ok){throw new RuntimeException("I/O error reading "+path);}}
		if(!header || out.size()!=expectedFamilies){throw new IllegalArgumentException("Roster "+path+" has "+out.size()+" families; expected "+expectedFamilies);}
		return out;
	}

	private enum Outcome {
		ASSIGNED(1,"ASSIGNED"), BAIT_INTERCEPTED(2,"BAIT_INTERCEPTED"),
		MASKED_INTERCEPTED(3,"MASKED_INTERCEPTED"), NO_PRESELECT_SURVIVOR(4,"NO_PRESELECT_SURVIVOR"),
		NO_ACCEPTANCE_PASS(5,"NO_ACCEPTANCE_PASS"), NO_VALID_KMER(6,"NO_VALID_KMER"),
		INVALID_SEQUENCE(7,"INVALID_SEQUENCE");
		final int code; final String text;
		Outcome(int code_,String text_){code=code_;text=text_;}
		static Outcome parse(LineParser1 lp,int term,String path,long lineNo){
			for(Outcome o : values()){if(lp.termEquals(o.text,term)){return o;}}
			throw bad(path,lineNo,"unknown assignment outcome");
		}
		static Outcome fromEncoded(int encoded){
			final int code=encoded>>>16;
			for(Outcome o : values()){if(o.code==code){return o;}}
			throw new IllegalArgumentException("Unknown encoded assignment outcome: "+encoded);
		}
	}

	private static AssignmentIndex loadAssignments(final String manifestPath, final HashRankMap roster){
		final AssignmentManifest manifest=AssignmentManifest.load(manifestPath);
		if(manifest.queries>Integer.MAX_VALUE){throw new IllegalArgumentException("Query count exceeds int map capacity: "+manifest.queries);}
		final AssignmentIndex out=new AssignmentIndex((int)manifest.queries);
		for(AssignmentFile file : manifest.files){
			if(!ReferenceCdsSurvivalLabelReader.sha256File(file.queryPath).equals(file.querySha256)){throw new IllegalArgumentException("Query-file SHA-256 mismatch: "+file.queryPath);}
			if(!sha80(file.hitPath).equals(file.hitSha80)){throw new IllegalArgumentException("Hit-file sha80 mismatch: "+file.hitPath);}
			if(!sha80(file.provenancePath).equals(file.provenanceSha80)){throw new IllegalArgumentException("Provenance-file sha80 mismatch: "+file.provenancePath);}
			final long before=out.size();
			final ProvenanceCounts counts=readProvenance(file.provenancePath,roster,out);
			if(counts.queries!=file.queries || counts.assigned!=file.assigned){throw new IllegalArgumentException("Assignment file "+file.provenancePath+" observed queries/assigned="+counts.queries+"/"+counts.assigned+" but manifest declares "+file.queries+"/"+file.assigned);}
			if(out.size()-before!=file.queries){throw new AssertionError("Assignment map grew by "+(out.size()-before)+" for "+file.queries+" queries");}
			final long hitRows=validateHitFile(file.hitPath,roster,out);
			if(hitRows!=file.assigned){throw new IllegalArgumentException("Hit file "+file.hitPath+" has "+hitRows+" rows; manifest declares "+file.assigned);}
			out.hitRows+=hitRows;
		}
		if(out.size()!=manifest.queries || out.hitRows!=manifest.assigned){throw new AssertionError("Loaded assignments/hits "+out.size()+"/"+out.hitRows+" != manifest "+manifest.queries+"/"+manifest.assigned);}
		return out;
	}

	private static long validateHitFile(final String path, final HashRankMap roster,
			final AssignmentIndex assignments){
		final ByteFile input=ByteFile.makeByteFile(path,true); final LineParser1 lp=new LineParser1((byte)'\t');
		long rows=0,lineNo=0; boolean ok=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lineNo++; if(line.length==0 || line[0]=='#'){continue;} lp.set(line);
				if(lp.terms()!=3 || !lp.termEquals("1.0",2)){throw bad(path,lineNo,"gene-hit row must be query, rep_id, literal 1.0");}
				lp.setBounds(0); final int encoded=assignments.peek(line,lp.a(),lp.b());
				if(encoded<0 || Outcome.fromEncoded(encoded)!=Outcome.ASSIGNED){throw bad(path,lineNo,"gene-hit query is absent or not ASSIGNED in provenance");}
				final int expectedRank=(encoded&0xFFFF)-1;
				lp.setBounds(1); final int observedRank=roster.get(line,lp.a(),lp.b());
				if(observedRank!=expectedRank){throw bad(path,lineNo,"gene-hit rep_id disagrees with provenance family_idx");}
				rows++;
			}
			ok=true;
		}finally{if(input.close() && ok){throw new RuntimeException("I/O error reading "+path);}}
		return rows;
	}

	private static ProvenanceCounts readProvenance(final String path, final HashRankMap roster, final AssignmentIndex out){
		final ByteFile input=ByteFile.makeByteFile(path,true); final LineParser1 lp=new LineParser1((byte)'\t');
		final ProvenanceCounts counts=new ProvenanceCounts();
		long lineNo=0; boolean header=false,ok=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lineNo++; if(line.length==0){continue;} lp.set(line);
				if(!header){
					if(lp.terms()!=11 || !lp.termEquals("#gene",0) || !lp.termEquals("assigned",1) || !lp.termEquals("rep_id",2) || !lp.termEquals("family_idx",3) || !lp.termEquals("reason",10)){throw bad(path,lineNo,"invalid assignment-provenance header");}
					header=true; continue;
				}
				if(lp.terms()!=11){throw bad(path,lineNo,"assignment-provenance row must have 11 fields");}
				final Outcome outcome=Outcome.parse(lp,10,path,lineNo);
				final boolean assigned=lp.termEquals("t",1);
				if(assigned!=(outcome==Outcome.ASSIGNED)){throw bad(path,lineNo,"assigned flag disagrees with outcome "+outcome.text);}
				int rank=-1;
				if(outcome==Outcome.ASSIGNED || outcome==Outcome.MASKED_INTERCEPTED){
					if(lp.termEquals("NA",2) || lp.termEquals("NA",3)){throw bad(path,lineNo,outcome.text+" lacks rep_id/family_idx");}
					rank=lp.parseInt(3);
					lp.setBounds(2); final int rosterRank=roster.get(line,lp.a(),lp.b());
					if(rank<0 || rank!=rosterRank){throw bad(path,lineNo,"rep_id/family_idx do not match the active roster");}
				}else if(!lp.termEquals("NA",3)){
					throw bad(path,lineNo,"non-assigned/non-masked outcome carries a family_idx");
				}
				if(rank>=0xFFFF){throw bad(path,lineNo,"family_idx exceeds contribution encoding");}
				final int encoded=(outcome.code<<16)+(outcome==Outcome.ASSIGNED ? rank+1 : 0);
				lp.setBounds(0); out.put(line,lp.a(),lp.b(),encoded,path+":"+lineNo);
				counts.queries++; if(outcome==Outcome.ASSIGNED){counts.assigned++;}
			}
			ok=true;
		}finally{if(input.close() && ok){throw new RuntimeException("I/O error reading "+path);}}
		if(!header){throw new IllegalArgumentException("Missing assignment-provenance header: "+path);}
		return counts;
	}

	private static final class ProvenanceCounts { long queries,assigned; }

	private static final class AssignmentFile {
		final int rank; final String queryPath,querySha256,hitPath,hitSha80,provenancePath,provenanceSha80; final long assigned,queries;
		AssignmentFile(int rank_,String queryPath_,String querySha256_,String hitPath_,String hitSha80_,long assigned_,String provenancePath_,String provenanceSha80_,long queries_){
			rank=rank_;queryPath=queryPath_;querySha256=querySha256_;hitPath=hitPath_;hitSha80=hitSha80_;assigned=assigned_;provenancePath=provenancePath_;provenanceSha80=provenanceSha80_;queries=queries_;
		}
	}

	private static final class AssignmentManifest {
		final ArrayList<AssignmentFile> files=new ArrayList<AssignmentFile>(); long assigned=0,queries=0;
		static AssignmentManifest load(final String path){
			final AssignmentManifest out=new AssignmentManifest();
			final File base=new File(path).getAbsoluteFile().getParentFile();
			final ByteFile input=ByteFile.makeByteFile(path,true);
			final LineParser1 lp=new LineParser1((byte)'\t');
			long lineNo=0, endChunks=-1, endAssigned=-1, endQueries=-1; boolean schema=false,columns=false,end=false,ok=false;
			try{
				for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
					lineNo++;
					if(line.length==0){continue;}
					if(end){throw bad(path,lineNo,"content after #end");}
					lp.set(line);
					if(line[0]=='#'){
						if(lp.termEquals("#schema_version",0)){if(schema || lp.terms()!=2 || !lp.termEquals(ASSIGNMENT_MANIFEST_SCHEMA,1)){throw bad(path,lineNo,"invalid schema");}schema=true;}
						else if(lp.termEquals("#columns",0)){if(columns || lp.terms()!=10 || !lp.termEquals("chunk_rank",1) || !lp.termEquals("query_path",2) || !lp.termEquals("query_sha256",3) || !lp.termEquals("hit_path",4) || !lp.termEquals("hit_sha80",5) || !lp.termEquals("assigned",6) || !lp.termEquals("provenance_path",7) || !lp.termEquals("provenance_sha80",8) || !lp.termEquals("queries",9)){throw bad(path,lineNo,"invalid columns");}columns=true;}
						else if(lp.termEquals("#end",0)){
							if(lp.terms()!=4 || !lp.termStartsWith("chunks=",1) || !lp.termStartsWith("assigned=",2) || !lp.termStartsWith("queries=",3)){throw bad(path,lineNo,"invalid trailer");}
							endChunks=parseLong(lp,1,7,path,lineNo); endAssigned=parseLong(lp,2,9,path,lineNo); endQueries=parseLong(lp,3,8,path,lineNo); end=true;
						}
						continue;
					}
					if(!schema || !columns || lp.terms()!=9){throw bad(path,lineNo,"manifest row before schema/columns or wrong field count");}
					final int rank=lp.parseInt(0);
					if(rank!=out.files.size()){throw bad(path,lineNo,"chunk_rank "+rank+" != contiguous row "+out.files.size());}
					final String rawQuery=lp.parseString(1), queryPath=new File(rawQuery).isAbsolute() ? rawQuery : new File(base,rawQuery).getPath();
					final String queryHash=lp.parseString(2); if(!ReferenceCdsShredSurvivalTableReader.isHex64(queryHash)){throw bad(path,lineNo,"invalid query_sha256");}
					final String rawHit=lp.parseString(3), hitPath=new File(rawHit).isAbsolute() ? rawHit : new File(base,rawHit).getPath();
					final String hitHash=lp.parseString(4); if(!isSha80(hitHash)){throw bad(path,lineNo,"invalid hit_sha80");}
					final long assigned=parseLong(lp,5,0,path,lineNo);
					final String rawProvenance=lp.parseString(6), provenancePath=new File(rawProvenance).isAbsolute() ? rawProvenance : new File(base,rawProvenance).getPath();
					final String provenanceHash=lp.parseString(7); if(!isSha80(provenanceHash)){throw bad(path,lineNo,"invalid provenance_sha80");}
					final long queries=parseLong(lp,8,0,path,lineNo); if(assigned>queries){throw bad(path,lineNo,"assigned exceeds queries");}
					out.files.add(new AssignmentFile(rank,queryPath,queryHash,hitPath,hitHash,assigned,provenancePath,provenanceHash,queries));
					out.assigned+=assigned; out.queries+=queries;
				}
				ok=true;
			}finally{if(input.close() && ok){throw new RuntimeException("I/O error reading "+path);}}
			if(!schema || !columns || !end || endChunks!=out.files.size() || endAssigned!=out.assigned || endQueries!=out.queries){throw new IllegalArgumentException("Incomplete or inconsistent assignment manifest: "+path);}
			return out;
		}
	}

	private static long parseLong(final LineParser1 lp, final int term, final int offset,
			final String path, final long lineNo){
		lp.setBounds(term); final byte[] line=lp.line(); final int from=lp.a()+offset, to=lp.b();
		if(from>=to){throw bad(path,lineNo,"empty integer term "+term);}
		long value=0;
		for(int i=from; i<to; i++){
			final byte c=line[i]; if(c<'0' || c>'9'){throw bad(path,lineNo,"non-digit in integer term "+term);}
			if(value>(Long.MAX_VALUE-9)/10){throw bad(path,lineNo,"integer overflow in term "+term);}
			value=value*10+(c-'0');
		}
		return value;
	}

	private static IllegalArgumentException bad(String path,long line,String message){return new IllegalArgumentException(path+":"+line+": "+message);}

	private static boolean isSha80(final String s){
		if(s==null || s.length()!=20){return false;}
		for(int i=0; i<s.length(); i++){final char c=s.charAt(i); if(!((c>='0' && c<='9') || (c>='a' && c<='f'))){return false;}}
		return true;
	}

	static String sha80(final String path){
		final String full=ReferenceCdsSurvivalLabelReader.sha256File(path);
		return full.substring(full.length()-20);
	}

	private static long hash1(final byte[] data, final int from, final int to){
		long h=0xcbf29ce484222325L;
		for(int i=from; i<to; i++){h^=(data[i]&0xFF); h*=0x100000001b3L;}
		return h;
	}

	private static long hash2(final byte[] data, final int from, final int to){
		long h=0x9e3779b97f4a7c15L^(to-from);
		for(int i=from; i<to; i++){
			h^=(data[i]&0xFF)+0x9e3779b97f4a7c15L;
			h=Long.rotateLeft(h,27)*0x94d049bb133111ebL;
		}
		return h^(h>>>31);
	}

	/** Avoid a magic narrowing check against {@link Integer#MAX_VALUE}. */
	private static final class TotalInt { static final long MAX=Integer.MAX_VALUE; }
}
