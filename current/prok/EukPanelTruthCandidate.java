package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Materializes policy v1.4 as a review candidate from immutable descriptive loci.
 * No locus reconstruction, truth freeze, caller execution or biological relabeling.
 * Existing members and source-only masks remain byte-preserved producer companions.
 * @author Raiden */
public final class EukPanelTruthCandidate {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, EukPanelTruthCandidate.class, false);
		String loci=null, members=null, links=null, sources=null, models=null, cases=null, diagnostics=null, out=null;
		int nl=-1, nm=-1, na=-1, nc=-1, nu=-1;long nh=-1, ne=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("loci")){loci=v;}else if(k.equals("members")){members=v;}
			else if(k.equals("links")){links=v;}else if(k.equals("sources")){sources=v;}
			else if(k.equals("models")){models=v;}else if(k.equals("cases")){cases=v;}
			else if(k.equals("diagnostics")){diagnostics=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("expectedloci")){nl=Integer.parseInt(v);}else if(k.equals("expectedmodels")){nm=Integer.parseInt(v);}
			else if(k.equals("expectedmembers")){nh=Long.parseLong(v);}else if(k.equals("expectedlinks")){ne=Long.parseLong(v);}
			else if(k.equals("expectedambiguous")){na=Integer.parseInt(v);}else if(k.equals("expectedcrossfold")){nc=Integer.parseInt(v);}
			else if(k.equals("expectedunexpected")){nu=Integer.parseInt(v);}else{throw new IllegalArgumentException("Unknown flag: "+k);}
		}
		require(loci!=null && members!=null && links!=null && sources!=null && models!=null && cases!=null && diagnostics!=null && out!=null
			&& nl>=0 && nm>0 && nh>=0 && ne>=0 && na>=0 && nc>=0 && nu>=0, "Explicit frozen inputs and denominators required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		final EukPanelTruthCandidate x=new EukPanelTruthCandidate();
		x.validation.loadLoci(loci);x.validation.loadMembers(members, nh);x.loadRows(loci, nl);
		x.loadModels(models, nm);x.loadHistory(links, EukPanelGffValidation.hashes(sources), ne);
		x.loadCases(cases, na);x.loadDiagnostics(diagnostics, nc, nu);x.write(out);
		Shared.closeStream(pp.outstream);
	}
	void loadRows(String path, int expected){
		assert(rows.isEmpty()):"One accepted descriptive generation supplies candidate identity and geometry";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in, EukPanelLoci.LOCUS_HEADER, path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==26, "Expected descriptive locus schema");
			final EukPanelLocusDiagnostics.Locus locus=validation.byId.get(p.parseString(0));
			require(locus!=null, "Candidate must retain a validated locus");final Row row=new Row(line, p, locus);
			require(byId.put(locus.id, row)==null, "Duplicate candidate locus");rows.add(row);
			for(EukPanelLoci.Hit h:locus.members){require(byHit.put(hitKey(h), h)==null, "Duplicate original hit identity");}
			require(seedRows.put(hitKey(locus.seed), row)==null, "One seed belongs to exactly one locus");
		}}finally{require(!in.close(), "Locus metadata read failed");}
		require(rows.size()==expected && rows.size()==validation.loci.size(), "Candidate locus denominator changed");
	}
	void loadModels(String path, int expected){
		assert(modelFamilies.isEmpty()):"The pinned model table, not naming guesses, defines target folds";
		final EukPanelLoci roster=new EukPanelLoci();roster.loadModels(path, expected);
		for(EukPanelLoci.Model m:roster.models.values()){if(m!=null){
			require(NcrnaClippingComparison.known(m.family) && modelFamilies.put(m.accession, m.family)==null, "Unknown or duplicate target accession");
		}}
		final LineParser1 p=new LineParser1('\t');
		for(EukPanelLoci.Hit h:byHit.values()){if(h!=null){
			final EukPanelLoci.Model m=roster.models.get(h.h.bundle+'\t'+h.h.model);
			p.set(h.h.raw);
			require(m!=null && m.family.equals(h.h.family) && m.clen==h.h.clen && m.ga==p.parseDouble(26), "Member disagrees with pinned model map");
		}}
	}
	void loadHistory(String path, ObjectMap<String,String> sources, long expected){
		assert(historyLinks==0):"Source links are physical evidence only; one pass must conserve the frozen link count";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in, LINK_HEADER, path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==38 && p.termEquals("PHYSICAL_LINK_REVIEW_ONLY",37), "Expected original physical source-link schema");
			final String assembly=p.parseString(1), run=p.parseString(27), source=p.parseString(2);
			final long hitLine=p.parseLong(28), sourceLine=p.parseLong(4);
			final String key=assembly+'\t'+run+'\t'+hitLine;final EukPanelLoci.Hit h=byHit.get(key);
			require(h!=null && p.termEquals("VALID_INTERVAL",5) && p.termEquals(h.h.nativeSeqid,8) && p.termEquals(h.h.strand,20)
				&& p.termEquals(h.h.model,29) && p.termEquals(h.h.family,30) && p.termEquals(h.h.strand,33), "Historical link must retain member identity and strand");
			require(p.parseString(3).equals(sources.get(source)) && sourceLine>0 && seen.add(source+'\t'+sourceLine+'\t'+key), "Unbound or repeated historical source link");
			final long start=p.parseLong(17), end=p.parseLong(18), from=p.parseLong(31), to=p.parseLong(32);
			require(start>0 && end>=start && Math.min(from,to)==h.h.start && Math.max(from,to)==h.h.stop, "Historical link coordinates changed");
			final long overlap=EukPanelCmEvidence.overlap(start,end,h.h.start,h.h.stop);
			require(overlap>0 && overlap==p.parseLong(34), "Historical support needs positive physical overlap");
			h.h.sourceLinks++;historyLinks++;
			// Seed-only links are exactly the annotations overlapping the fixed locus span.
			// A member can extend beyond that span, so its annotation alone is insufficient.
			final Row row=seedRows.get(key);final String type=p.parseString(16);
			if(row!=null && type.endsWith("RNA") && !type.equals("mRNA")){
				row.histAny=true;final String family=sourceFamily(type,p.parseString(22));
				if(row.locus.family.equals(family)){row.histSame=true;}
			}
		}}finally{require(!in.close(), "Historical link read failed");}
		require(historyLinks==expected, "Historical link denominator changed");
		for(EukPanelLoci.Hit h:byHit.values()){if(h!=null){require(h.h.sourceLinks==h.h.expectedLinks, "Historical member link conservation failed: "+hitKey(h));}}
	}
	public String sourceFamily(String type, String attrs){
		assert(type.endsWith("RNA") && !type.equals("mRNA")):"Historical support excludes mRNA, genes, regions and exons";
		String family=null;final String db=EukPanelCmProposal.gffAttribute(attrs,"Dbxref");
		if(db!=null){for(String token:db.split(",")){if(token.startsWith("RFAM:")){
			final String f=modelFamilies.get(token.substring(5));if(f==null || family!=null && !family.equals(f)){return "UNKNOWN";}family=f;
		}}}
		if(type.equals("tRNA") || type.equals("tmRNA")){return family==null || family.equals(type)?type:"UNKNOWN";}
		if(family!=null){return family;}
		if(!type.equals("rRNA")){return "UNKNOWN";}
		final String product=EukPanelCmProposal.gffAttribute(attrs,"product");if(product==null){return "UNKNOWN";}
		final String text=product.toLowerCase(Locale.ROOT);
		final String[] names={"5.8s","5s","16s","18s","23s","25s","26s","28s"};
		final String[] folds={"5.8S","5S","SSU","SSU","LSU","LSU","LSU","LSU"};
		for(int i=0;i<names.length;i++){if(word(text,names[i])){if(family!=null && !family.equals(folds[i])){return "UNKNOWN";}family=folds[i];}}
		return family==null?"UNKNOWN":family;
	}
	static boolean word(String text, String token){
		assert(!token.isEmpty()):"Explicit rRNA names require token boundaries, not length-based family guesses";
		for(int at=text.indexOf(token);at>=0;at=text.indexOf(token,at+1)){
			final int end=at+token.length();if((at==0 || !Character.isLetterOrDigit(text.charAt(at-1)) && text.charAt(at-1)!='.')
				&& (end==text.length() || !Character.isLetterOrDigit(text.charAt(end)))){return true;}
		}return false;
	}
	void loadCases(String path, int expected){
		assert(ambiguous.isEmpty()):"Every spanning digest binds the existing locus and two original witnesses";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in, CASE_HEADER, path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==18, "Spanning case schema changed");final Row row=byId.get(p.parseString(0));
			require(row!=null && row.ambiguous && p.termEquals(row.locus.assembly,2) && p.termEquals(row.locus.seqid,3) && p.termEquals(row.locus.strand,4)
				&& p.termEquals(row.locus.family,5) && p.termEquals(row.locus.run,6) && p.parseLong(7)==row.locus.line
				&& p.parseLong(8)==row.locus.start && p.parseLong(9)==row.locus.end, "Spanning case differs from its locus");
			final EukPanelLoci.Hit a=witness(row,p,10), b=witness(row,p,14);
			require(EukPanelLoci.substantial(a) && EukPanelLoci.substantial(b) && EukPanelLoci.separate(a,b)
				&& EukPanelLoci.overlap(row.locus.seed,a)>0 && EukPanelLoci.overlap(row.locus.seed,b)>0, "Spanning witnesses must be separate GA/full hits overlapping the seed");
			require(p.termEquals("SPANNING_SEED",1) || p.termEquals("SEPARATE_MEMBERS",1), "Unknown spanning reason");
			if(ambiguous.add(row.locus.id)){row.witnessA=a;row.witnessB=b;}
		}}finally{require(!in.close(), "Spanning cases read failed");}
		int count=0;for(Row row:rows){if(row.ambiguous){require(ambiguous.contains(row.locus.id), "Flagged locus lacks a witness digest");count++;}}
		require(count==expected && ambiguous.size()==expected, "Ambiguous locus denominator changed");
	}
	EukPanelLoci.Hit witness(Row row, LineParser1 p, int column){
		assert(column==10 || column==14):"Witness offsets follow the immutable 18-column case table";
		final EukPanelLoci.Hit h=byHit.get(row.locus.assembly+'\t'+p.parseString(column)+'\t'+p.parseLong(column+1));
		require(h!=null && h.partition.equals(row.locus.seed.partition) && h.h.start==p.parseLong(column+2) && h.h.stop==p.parseLong(column+3), "Spanning witness identity/coordinates changed");return h;
	}
	void loadDiagnostics(String path, int expectedCross, int expectedUnexpected){
		assert(crossfold.isEmpty() && unexpected.isEmpty()):"Diagnostic cases are consumed once; flags cannot silently veto admission";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in, DIAGNOSTIC_HEADER, path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==35 && p.termEquals("REVIEW_ONLY",34), "Diagnostic schema changed");
			final Row a=byId.get(p.parseString(1)), b=p.termEquals("NA",2)?null:byId.get(p.parseString(2));
			require(a!=null && a.locus.ga, "Diagnostic needs GA-supported locus");
			final String reason=p.parseString(0);require(seen.add(reason+'\t'+p.parseString(1)+'\t'+p.parseString(2)), "Repeated diagnostic case");
			require(p.termEquals(a.locus.assembly,3) && p.termEquals(a.locus.seqid,4) && p.termEquals(a.locus.strand,5) && p.termEquals(a.locus.compartment,6), "Diagnostic locus identity changed");
			final EukPanelLoci.Hit h=diagnosticHit(a,p,8);
			if(reason.equals("AMBIGUOUS_CROSSFOLD")){
				require(b!=null && b.locus.ga && EukPanelLocusDiagnostics.ribosomal(a.locus) && EukPanelLocusDiagnostics.ribosomal(b.locus)
					&& !a.locus.family.equals(b.locus.family) && a.locus.partition.equals(b.locus.partition), "Crossfold pair must be two same-strand supported ribosomal folds");
				final EukPanelLoci.Hit z=diagnosticHit(b,p,21);
				require(h==a.locus.seed && z==b.locus.seed && EukPanelLoci.overlap(h,z)>0 && EukPanelLoci.overlap(h,z)==p.parseLong(7), "Crossfold case must use fixed seed spans");
				a.crossfold=b.crossfold=true;crossfold.add(new Pair(a,b,h,z));
			}else if(reason.equals("COMPARTMENT_UNEXPECTED")){
				require(b==null && a.locus.compartment.equals("NUCLEAR") && h.h.gaPass && (h.h.model.equals("RF00010") || h.h.model.equals("RF00023")), "Unexpected compartment definition changed");
				require(unexpected.put(a.locus.id,h)==null, "Expected one supporting-hit digest per unexpected locus");a.unexpected=true;
			}else if(reason.equals("SEED_BELOW_GA")){
				require(b==a && !a.locus.seedGa && h==a.locus.seed && diagnosticHit(a,p,21).h.gaPass, "Below-GA seed needs a supported member");
			}else{throw new IllegalArgumentException("Unknown diagnostic reason: "+reason);}
		}}finally{require(!in.close(), "Diagnostic read failed");}
		require(crossfold.size()==expectedCross && unexpected.size()==expectedUnexpected, "Crossfold/unexpected denominator changed");
	}
	EukPanelLoci.Hit diagnosticHit(Row row, LineParser1 p, int c){
		assert(c==8 || c==21):"Each diagnostic side has thirteen fields in the validated case schema";
		final EukPanelLoci.Hit h=byHit.get(row.locus.assembly+'\t'+p.parseString(c+4)+'\t'+p.parseLong(c+6));
		require(h!=null && row.locus.members.contains(h) && p.termEquals(h.h.family,c) && p.termEquals(h.h.model,c+1)
			&& p.parseLong(c+2)==h.h.start && p.parseLong(c+3)==h.h.stop && p.termEquals(h.table,c+5) && p.parseDouble(c+7)==h.bits
			&& p.parseInt(c+9)==h.h.mf && p.parseInt(c+10)==h.h.mt && p.parseInt(c+11)==h.h.clen && p.termEquals(h.h.trunc,c+12), "Diagnostic witness differs from original member");return h;
	}
	void write(String prefix){
		assert(rows.size()==validation.loci.size()):"Candidate output cannot remove a descriptive locus";
		RrnaResourceIO.fresh(prefix,".truth_candidate.tsv.gz",".candidate_counts.tsv",".spanning_digest.tsv",".crossfold_digest.tsv",".unexpected_digest.tsv");
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".truth_candidate.tsv.gz"), counts=RrnaResourceIO.writer(prefix+".candidate_counts.tsv");
		final ByteStreamWriter spanning=RrnaResourceIO.writer(prefix+".spanning_digest.tsv"), cross=RrnaResourceIO.writer(prefix+".crossfold_digest.tsv"), odd=RrnaResourceIO.writer(prefix+".unexpected_digest.tsv");
		final ByteBuilder b=new ByteBuilder();final ObjectMap<String,long[]> bins=new ObjectMap<String,long[]>(String.class,long[].class);long members=0;
		try{
			out.println(CANDIDATE_HEADER);spanning.println(DIGEST_HEADER);cross.println(DIGEST_HEADER);odd.println(DIGEST_HEADER);
			for(Row row:rows){
				final EukPanelLocusDiagnostics.Locus x=row.locus;final String admission=x.ga?"ADMITTED":"BELOW_GA";
				// All 25 declared model accessions belong to shipped folds. Source-only
				// OUT_OF_MAP entries stay in their separate provenance-bound mask.
				final String flags=row.flags();out.print(b.clear().append(row.raw).tab().append(admission).tab().append(row.geometry).tab().append(flags)
					.tab().append(row.histAny?1:0).tab().append(row.histSame?1:0).tab().append(x.ga && x.seedGa?1:0).tab().append("FREEZE_CANDIDATE_NOT_FROZEN").nl());
				EukPanelCmProposal.increment(bins,row.bin+'\t'+x.family+'\t'+x.compartment+'\t'+admission+'\t'+row.geometry,1);members+=x.members.size();
				if(row.ambiguous){digest(spanning,b,"SPANNING_AMBIGUOUS",row,null,x.seed,row.witnessA,row.witnessB,"REVIEW_SEPARATE_WITNESSES_KEEP_FLAGGED_CANDIDATE");}
				if(row.unexpected){digest(odd,b,"COMPARTMENT_UNEXPECTED",row,null,unexpected.get(x.id),null,null,"RETAIN_ADMITTED_FLAGGED_NO_COMPARTMENT_RECLASSIFICATION");}
			}
			for(Pair pair:crossfold){digest(cross,b,"AMBIGUOUS_CROSSFOLD",pair.a,pair.b,pair.h,pair.z,null,"RETAIN_BOTH_ADMITTED_FLAGGED_PENDING_OWNER_REVIEW");}
			counts.println("selection_bin\tfold_family\tactual_compartment\tadmission_class\tgeometry_stratum\tloci");long total=0;
			for(String key:EukPanelCmProposal.sorted(bins)){final long n=bins.get(key)[0];counts.print(b.clear().append(key).tab().append(n).nl());total+=n;}
			require(total==rows.size() && members==validation.memberCount, "Candidate bins and members must conserve the validated denominators");
		}finally{RrnaResourceIO.finish(out,counts,spanning,cross,odd);}
		System.err.println("TRUTH_CANDIDATE_COMPLETE loci="+rows.size()+" members="+members+" links="+historyLinks+" spanning="+ambiguous.size()+" crossfold="+crossfold.size()+" unexpected="+unexpected.size()+" frozen=false");
	}
	static void digest(ByteStreamWriter out, ByteBuilder b, String reason, Row a, Row z, EukPanelLoci.Hit h, EukPanelLoci.Hit w1, EukPanelLoci.Hit w2, String recommendation){
		assert(h!=null):"Every readable case carries a real original hit witness";
		b.clear().append(reason).tab().append(a.locus.id).tab().append(z==null?"NA":z.locus.id).tab().append(a.locus.assembly).tab().append(a.locus.seqid)
			.tab().append(a.locus.strand).tab().append(a.locus.compartment).tab();hitDigest(b,h);b.tab();hitDigest(b,w1);b.tab();hitDigest(b,w2);
		out.print(b.tab().append(recommendation).tab().append("RECOMMENDATION_ONLY").nl());
	}
	static void hitDigest(ByteBuilder b, EukPanelLoci.Hit h){
		assert(b!=null):"One reusable byte buffer carries a bounded case line";
		if(h==null){b.append("NA");return;}
		b.append(h.h.model).append(':').append(h.h.start).append('-').append(h.h.stop).append(";bits=").append(h.bitsText)
			.append(";model=").append(h.h.mf).append('-').append(h.h.mt).append('/').append(h.h.clen).append(";GA=").append(h.h.gaPass?1:0)
			.append(";trunc=").append(h.h.trunc).append(";source=").append(h.h.run).append(':').append(h.h.lineNumber);
	}
	static String hitKey(EukPanelLoci.Hit h){return h.h.assembly+'\t'+h.h.run+'\t'+h.h.lineNumber;}
	public static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Row{
		Row(byte[] line, LineParser1 p, EukPanelLocusDiagnostics.Locus x){
			raw=line;locus=x;bin=p.parseString(7);ambiguous=p.termEquals("SPANNING_AMBIGUOUS",22);
			require(ambiguous || p.termEquals("NONE",22), "Unknown descriptive flag");
			final boolean full=EukPanelLocusDiagnostics.bit(p,15), truncated=EukPanelLocusDiagnostics.bit(p,16);final String edge=p.parseString(17);
			require(full==x.seed.h.full && truncated==!x.seed.h.trunc.equals("no") && edge.equals(x.seed.h.edge), "Seed geometry columns disagree with original member");
			geometry=truncated || edge.equals("LEFT") || edge.equals("RIGHT") || edge.equals("BOTH")?"EDGE":edge.equals("UNASSESSABLE")?"UNASSESSABLE":full?"FULLSPAN_INTACT":"PARTIAL_INTERNAL";
			final LineParser1 hit=new LineParser1('\t');hit.set(x.seed.h.raw);crossFamily=hit.parseLong(38)>0;
		}
		String flags(){
			assert(!histSame || histAny):"Same-family historical RNA support implies any-RNA support";
			final ByteBuilder b=new ByteBuilder();if(ambiguous){flag(b,"SPANNING_AMBIGUOUS");}if(crossFamily){flag(b,"CROSS_FAMILY_OVERLAP");}
			if(crossfold){flag(b,"AMBIGUOUS_CROSSFOLD");}if(unexpected){flag(b,"COMPARTMENT_UNEXPECTED");}if(locus.ga && !locus.seedGa){flag(b,"SEED_BELOW_GA");}
			flag(b,histAny?"HIST_SUPPORTED":"HIST_UNSUPPORTED");if(histAny){flag(b,"HIST_ANY_RNA");}if(histSame){flag(b,"HIST_SAME_FAMILY");}return b.toString();
		}
		static void flag(ByteBuilder b,String flag){if(b.length()>0){b.append(';');}b.append(flag);}
		final byte[] raw;final EukPanelLocusDiagnostics.Locus locus;final String bin,geometry;final boolean ambiguous,crossFamily;
		boolean crossfold,unexpected,histAny,histSame;EukPanelLoci.Hit witnessA,witnessB;
	}
	static final class Pair{Pair(Row a_,Row b_,EukPanelLoci.Hit h_,EukPanelLoci.Hit z_){a=a_;b=b_;h=h_;z=z_;}final Row a,b;final EukPanelLoci.Hit h,z;}
	public static final String LINK_HEADER=EukPanelCmProposal.SOURCE_HEADER+"\tsource_srp_subtype\trun_id\ttable_line\tcm_model\tcm_family\tcm_seq_from\tcm_seq_to\tcm_strand\toverlap_bp\tsource_fraction\thit_fraction\tproposal_status";
	public static final String CASE_HEADER="locus_id\treason\tassembly\tnative_seqid\tstrand\tfold_family\tseed_run\tseed_line\tseed_start\tseed_end\twitness_a_run\twitness_a_line\twitness_a_start\twitness_a_end\twitness_b_run\twitness_b_line\twitness_b_start\twitness_b_end";
	static final String DIAGNOSTIC_HEADER="reason\tlocus_a\tlocus_b\tassembly\tnative_seqid\tstrand\tcompartment\toverlap_bp\tfamily_a\tmodel_a\tstart_a\tend_a\trun_a\ttable_a\tline_a\tbits_a\tga_a\tmodel_from_a\tmodel_to_a\tclen_a\ttrunc_a\tfamily_b\tmodel_b\tstart_b\tend_b\trun_b\ttable_b\tline_b\tbits_b\tga_b\tmodel_from_b\tmodel_to_b\tclen_b\ttrunc_b\tstatus";
	public static final String CANDIDATE_HEADER=EukPanelLoci.LOCUS_HEADER.replace("\tflags\tstatus\t","\tsource_flags\tsource_status\t")+"\tadmission_class\tgeometry_stratum\tflags\thist_any_rna\thist_same_family\tboundary_metric_eligible\tcandidate_status";
	static final String DIGEST_HEADER="reason\tlocus_a\tlocus_b\tassembly\tnative_seqid\tstrand\tcompartment\tseed_or_supporting_hit\twitness_a\twitness_b\trecommended_disposition\tstatus";
	final EukPanelLocusDiagnostics validation=new EukPanelLocusDiagnostics();final ArrayList<Row> rows=new ArrayList<Row>();
	final ObjectMap<String,Row> byId=new ObjectMap<String,Row>(String.class,Row.class),seedRows=new ObjectMap<String,Row>(String.class,Row.class);
	final ObjectMap<String,EukPanelLoci.Hit> byHit=new ObjectMap<String,EukPanelLoci.Hit>(String.class,EukPanelLoci.Hit.class),unexpected=new ObjectMap<String,EukPanelLoci.Hit>(String.class,EukPanelLoci.Hit.class);
	public final ObjectMap<String,String> modelFamilies=new ObjectMap<String,String>(String.class,String.class);
	final ObjectSet<String> ambiguous=new ObjectSet<String>(String.class);final ArrayList<Pair> crossfold=new ArrayList<Pair>();long historyLinks;
}
