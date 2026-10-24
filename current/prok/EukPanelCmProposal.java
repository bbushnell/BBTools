package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import parse.Parse;
import structures.ByteBuilder;

/** Descriptive, per-hit overlay for human truth review. Never admits or merges loci.
 * The producer verifies all consumed inputs against the accepted evidence seals.
 * Historical source and association tables remain immutable companion artifacts.
 * Exact model spans and contig-edge contact are observations, not partial-gene labels.
 * @author Ganyu */
public final class EukPanelCmProposal {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,EukPanelCmProposal.class,false);
		String inputs=null,sources=null,contigs=null,contigSources=null,hits=null,pairs=null,summary=null,features=null,associations=null,out=null,verified=null;
		boolean prokInputs=false;
		int expected=-1,expectedModels=-1,expectedRuns=-1;long expectedHits=-1,expectedSources=-1,expectedAssociations=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");
			final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("inputs")){inputs=v;}else if(k.equals("sources")){sources=v;}
			else if(k.equals("contigs")){contigs=v;}else if(k.equals("contigsources")){contigSources=v;}
			else if(k.equals("hits")){hits=v;}else if(k.equals("pairs")){pairs=v;}
			else if(k.equals("summary")){summary=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("features")){features=v;}else if(k.equals("associations")){associations=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}
			else if(k.equals("expectedmodels")){expectedModels=Integer.parseInt(v);}
			else if(k.equals("expectedruns")){expectedRuns=Integer.parseInt(v);}
			else if(k.equals("expectedhits")){expectedHits=Long.parseLong(v);}
			else if(k.equals("expectedsources")){expectedSources=Long.parseLong(v);}
			else if(k.equals("expectedassociations")){expectedAssociations=Long.parseLong(v);}
			else if(k.equals("prokinputs")){prokInputs=Parse.parseBoolean(v);}else if(k.equals("verified")){verified=v;}
			else{throw new IllegalArgumentException(k);}
		}
		require(inputs!=null && sources!=null && contigs!=null && contigSources!=null && hits!=null && pairs!=null && summary!=null && features!=null && associations!=null && out!=null
			&& expected>0 && expectedModels>0 && expectedRuns>0 && expectedHits>=0 && expectedSources>=0 && expectedAssociations>=0,"Explicit input bindings and expected denominators required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		final EukPanelCmProposal t=new EukPanelCmProposal();
		if(prokInputs){
			require(verified!=null,"Prok proposal requires producer-verified source bindings");
			final EukPanelCmEvidence binding=new EukPanelCmEvidence();binding.prokInputs=true;binding.verified=EukPanelGffValidation.hashes(verified);
			binding.panel.hashes=EukPanelGffValidation.hashes(sources);binding.panel.contigHashes=EukPanelGffValidation.hashes(contigSources);
			binding.loadProkInputs(inputs,expected);binding.loadProkContigs(contigs);t.panel=binding.panel;
			for(EukPanelGffValidation.Assembly a:t.panel.assemblies){require(t.bins.put(a.id,a.group)==null,"One prok reporting domain per assembly");}
		}else{
			require(verified==null,"verified= is scoped to explicit prokinputs=t");
			t.panel.hashes=EukPanelGffValidation.hashes(sources);t.panel.contigHashes=EukPanelGffValidation.hashes(contigSources);
			t.panel.loadInputs(inputs,expected);t.panel.loadContigs(contigs);t.loadBins(inputs);
		}
		t.loadHits(hits,expectedHits);t.loadPairs(pairs);t.loadSummary(summary,expectedModels,expectedRuns);
		t.write(out);t.writeSourceLinks(features,associations,out,expectedSources,expectedAssociations);Shared.closeStream(pp.outstream);
	}
	void loadBins(String path){
		assert(bins.isEmpty()):"One frozen selection supplies the reporting bins";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{p.set(in.nextLine());require(p.terms()==15 && p.termEquals("bin",2),"Expected selected-input schema");
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==15,"Input schema changed");
				final String id=p.parseString(1),bin=p.parseString(2);require(panel.byAssembly.get(id)!=null && !bin.isEmpty() && bins.put(id,bin)==null,"Unique known assembly/bin required");
			}
		}finally{require(!in.close(),"Input read failed");}
		require(bins.size()==panel.assemblies.size(),"Every selected assembly requires a reporting bin");
	}
	void loadHits(String path,long expected){
		assert(hits.isEmpty()):"Read each raw hit exactly once";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,HIT_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==29,"Expected normalized raw29 schema");
			final Hit h=new Hit(line,p);final EukPanelGffValidation.Assembly a=panel.byAssembly.get(h.assembly);
			require(a!=null && byHit.put(h.key,h)==null,"Unknown assembly or duplicate raw identity: "+h.key);
			h.bin=bins.get(h.assembly);
			final EukPanelGffValidation.Resolution r=EukPanelGffValidation.resolve(a,p.parseString(9));
			final String range=r.contig==null?"UNRESOLVED_SEQID":h.stop>r.contig.length?"OUT_OF_RANGE":"VALID_INTERVAL";
			require(p.termEquals(range,13) && p.termEquals(r.kind,12),"Hit mapping/range differs from accepted contig evidence");
			if(r.contig==null){require(h.nativeSeqid.equals("NA") && h.compartment.equals("UNKNOWN"),"Unresolved hits must stay explicit");h.edge="UNASSESSABLE";}
			else{
				require(h.nativeSeqid.equals(r.contig.id) && h.compartment.equals(r.contig.compartment),"Native identity/compartment mismatch");
				h.contigLength=r.contig.length;
				h.edge=!h.valid?"UNASSESSABLE":h.start==1?(h.stop==h.contigLength?"BOTH":"LEFT"):h.stop==h.contigLength?"RIGHT":"INTERNAL";
			}
			hits.add(h);increment(observed,h.summaryKey,1);
		}}finally{require(!in.close(),"Hit read failed");}
		require(hits.size()==expected,"Raw hit count differs from accepted evidence: "+hits.size()+"/"+expected);
	}
	/** Streaming lockstep join preserves source fields and exact physical associations.
	 * The all-feature source table is retained upstream. Here, all positive links
	 * are emitted, plus no-link non-mRNA RNA rows and all invalid source intervals.
	 * Neither region/gene/exon links nor contig contact become confirmed RNA loci. */
	void writeSourceLinks(String features,String associations,String prefix,long expectedSources,long expectedAssociations){
		assert(!scans.isEmpty()):"Source links can refer only to already validated scan hits";
		RrnaResourceIO.fresh(prefix,".source_links.tsv.gz",".source_review.tsv.gz");
		final ByteFile src=RrnaResourceIO.open(features),edges=RrnaResourceIO.open(associations);
		final ByteStreamWriter links=RrnaResourceIO.writer(prefix+".source_links.tsv.gz"),review=RrnaResourceIO.writer(prefix+".source_review.tsv.gz");
		final LineParser1 p=new LineParser1('\t'),q=new LineParser1('\t');final ByteBuilder b=new ByteBuilder();
		final ObjectMap<String,long[]> lastLines=new ObjectMap<String,long[]>(String.class,long[].class);
		long sourceRows=0,positive=0,unassociated=0;
		EukPanelGffValidation.Assembly current=null;long[] last=null;
		try{
			RrnaResourceIO.header(src,SOURCE_HEADER,features);RrnaResourceIO.header(edges,ASSOCIATION_HEADER,associations);
			links.println(SOURCE_HEADER+"\tsource_srp_subtype\trun_id\ttable_line\tcm_model\tcm_family\tcm_seq_from\tcm_seq_to\tcm_strand\toverlap_bp\tsource_fraction\thit_fraction\tproposal_status");
			review.println(SOURCE_HEADER+"\tsource_srp_subtype\tproposal_status");
			byte[] edge=edges.nextLine();
			for(byte[] line;(line=src.nextLine())!=null;){
				p.set(line);require(p.terms()==26,"Expected original26-column source evidence");
				if(current==null || !p.termEquals(current.id,1)){
					current=panel.byAssembly.get(p.parseString(1));require(current!=null,"Source assembly absent from selected panel");
					last=lastLines.get(current.id);if(last==null){last=new long[1];lastLines.put(current.id,last);}
				}
				final String assembly=current.id,gff=current.gff;final long sourceLine=p.parseLong(4);
				require(p.termEquals(gff,2) && p.termEquals(current.gffHash,3),"Source row must retain bound assembly/GFF identity");
				require(sourceLine>last[0],"Original source lines must remain unique and ordered within assembly");last[0]=sourceLine;
				final boolean valid=p.termEquals("VALID_INTERVAL",5),rna=nonMessengerRnaType(p);
				String srp=null;long matched=0,rows=0;
				while(edge!=null){
					q.set(edge);require(q.terms()==10,"Expected physical association schema");
					if(!q.termEquals(assembly,0) || !q.termEquals(gff,1) || q.parseLong(2)!=sourceLine){break;}
					rows++;
					if(q.termEquals("OVERLAP",3)){
						final Hit h=byHit.get(key(q.parseString(4),q.parseLong(5)));
						require(valid && h!=null && h.valid && h.assembly.equals(assembly) && p.termEquals(h.nativeSeqid,8) && p.termEquals(h.strand,20)
							&& q.termEquals(h.model,6),"Physical source/hit link identity disagrees");
						final long start=p.parseLong(17),stop=p.parseLong(18),n=EukPanelCmEvidence.overlap(start,stop,h.start,h.stop);
						require(start>0 && stop>=start && n>0 && n==q.parseLong(7),"Source/hit overlap must preserve both coordinate systems");
						require(Math.abs(q.parseDouble(8)-n/(double)(stop-start+1))<1e-8 && Math.abs(q.parseDouble(9)-n/(double)(h.stop-h.start+1))<1e-8,"Reciprocal physical-overlap fractions disagree");
						if(srp==null){srp=srpSubtype(p.parseString(22),rna);}
						b.clear().append(line).tab().append(srp).tab().append(h.run).tab().append(h.lineNumber).tab().append(h.model).tab().append(h.family)
							.tab().append(h.strand.equals("+")?h.start:h.stop).tab().append(h.strand.equals("+")?h.stop:h.start).tab().append(h.strand)
							.tab().append(n).tab().append(q.parseDouble(8),9).tab().append(q.parseDouble(9),9).tab().append("PHYSICAL_LINK_REVIEW_ONLY");links.print(b.nl());
						matched++;positive++;h.sourceLinks++;
					}else{require(sameTerm(p,23,q,3) && q.termEquals("NA",4),"Unassociated sentinel must match preserved source status");}
					edge=edges.nextLine();
				}
				require(rows>0 && matched==p.parseLong(24) && (matched>0?rows==matched:rows==1),"Each source must retain all physical links or one explicit sentinel");
				if((rna && matched==0) || !valid){
					if(srp==null){srp=srpSubtype(p.parseString(22),rna);}
					review.print(b.clear().append(line).tab().append(srp).tab().append("UNRESOLVED_SOURCE_REVIEW_ONLY").nl());unassociated++;
				}
				sourceRows++;
			}
			require(edge==null,"Association rows remain after source EOF");
			require(sourceRows==expectedSources && positive==expectedAssociations,"Source and physical-association totals differ from accepted evidence");
			for(Hit h:hits){require(h.sourceLinks==h.expectedLinks,"All raw-hit historical link counts must survive source join: "+h.key);}
		}finally{final boolean a=src.close(),z=edges.close();RrnaResourceIO.finish(links,review);require(!a && !z,"Source/association read failed");}
		System.err.println("CM_PROPOSAL_SOURCE_LINKS source_rows="+sourceRows+" positive_links="+positive+" unresolved_review_rows="+unassociated+" no_truth_admission=true");
	}
	static boolean nonMessengerRnaType(LineParser1 p){
		assert(p.terms()==26):"Source-type filtering uses the preserved source schema";
		if(p.termEquals("mRNA",16)){return false;}p.setBounds(16);final byte[] x=p.line();final int end=p.b();
		return end-p.a()>=3 && x[end-3]=='R' && x[end-2]=='N' && x[end-1]=='A';
	}
	static boolean sameTerm(LineParser1 a,int ai,LineParser1 b,int bi){
		assert(ai<a.terms() && bi<b.terms()):"Compare only declared source/association columns";
		a.setBounds(ai);b.setBounds(bi);final int n=a.b()-a.a();if(n!=b.b()-b.a()){return false;}
		for(int i=0;i<n;i++){if(a.line()[a.a()+i]!=b.line()[b.a()+i]){return false;}}return true;
	}
	public static String srpSubtype(String attrs,boolean rna){
		assert(attrs!=null):"Source attributes are preserved even when no subtype can be assigned";
		final String db=gffAttribute(attrs,"Dbxref");boolean small=false,large=false,other=false;
		if(db!=null){for(String token:db.split(",")){
			if(token.equals("RFAM:RF00169") || token.equals("RFAM:RF04183")){small=true;}
			else if(token.equals("RFAM:RF01854") || token.equals("RFAM:RF01857")){large=true;}
			else if(token.startsWith("RFAM:")){other=true;}
		}}
		// Review-only overlay: outside-map eukaryotic SRP and conflicting evidence
		// must survive intact (EUK_ORGANELLE_PANEL_TRUTH, SRP scope decision).
		// Never use product, length or contig geometry to choose small versus large.
		if((small && large) || (other && (small || large))){return "AMBIGUOUS_SRP_SUBTYPE";}
		if(small){return "SRP_small";}if(large){return "SRP_large";}
		final String product=gffAttribute(attrs,"product");
		final String lower=product==null?"":product.toLowerCase(Locale.ROOT);
		final String rnaClass=gffAttribute(attrs,"ncRNA_class");
		final boolean described=rna && (lower.contains("signal recognition particle") || lower.contains("srp") || "SRP_RNA".equals(rnaClass));
		if(described){return other?"OUTSIDE_CURRENT_SUBTYPE_MAP":"UNRESOLVED_SRP_SUBTYPE";}
		return "NOT_CLASSIFIED_AS_SRP";
	}
	/** GFF field boundaries follow GffLine.partial(): semicolons, never whitespace.
	 * Return original value bytes as text; callers never rewrite the attribute column. */
	public static String gffAttribute(String attrs,String key){
		assert(attrs!=null && key!=null && !key.isEmpty()):"Exact GFF key/value extraction requires intact attributes";
		final String prefix=key+'=';String value=null;
		for(int from=0,end=0;from<attrs.length();from=end+1){
			end=attrs.indexOf(';',from);if(end<0){end=attrs.length();}
			require(!(end-from==key.length() && attrs.startsWith(key,from)),"Malformed GFF attribute without value: "+key);
			if(end-from>=prefix.length() && attrs.startsWith(prefix,from)){
				require(value==null,"Duplicate GFF attribute: "+key);require(end>from+prefix.length(),"Malformed empty GFF attribute: "+key);
				value=attrs.substring(from+prefix.length(),end);
			}
		}
		return value;
	}
	void loadPairs(String path){
		assert(pairCount==0):"Every physical hit pair is counted once";
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,PAIR_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==11,"Expected physical-pair schema");
			final Hit a=byHit.get(key(p.parseString(3),p.parseLong(4))),b=byHit.get(key(p.parseString(6),p.parseLong(7)));
			require(a!=null && b!=null && a!=b && a.valid && b.valid,"Pair must name two valid retained hits");
			require(a.assembly.equals(b.assembly) && a.nativeSeqid.equals(b.nativeSeqid) && a.strand.equals(b.strand)
				&& p.termEquals(a.assembly,0) && p.termEquals(a.nativeSeqid,1) && p.termEquals(a.strand,2)
				&& p.termEquals(a.model,5) && p.termEquals(b.model,8),"Pair identities/strand/model must agree with raw hits");
			final long overlap=EukPanelCmEvidence.overlap(a.start,a.stop,b.start,b.stop);
			require(overlap>0 && overlap==p.parseLong(9),"Pair overlap must equal inclusive genomic overlap");
			final String pair=a.key.compareTo(b.key)<0?a.key+'\t'+b.key:b.key+'\t'+a.key;
			require(seen.add(pair),"Repeated/reversed pair inflates ambiguity counts");
			final boolean sameModel=a.model.equals(b.model);
			require(p.termEquals(sameModel?"SAME_MODEL":"CROSS_MODEL",10),"Pair relation disagrees with accession identity");
			if(sameModel){a.sameModel++;b.sameModel++;}else if(a.family.equals(b.family)){a.sameFamily++;b.sameFamily++;}
			else{a.crossFamily++;b.crossFamily++;}pairCount++;
		}}finally{require(!in.close(),"Pair read failed");}
	}
	void loadSummary(String path,int expectedModels,int expectedRuns){
		assert(scans.isEmpty()):"Keep zero-hit scans as a separate denominator";
		final ObjectSet<String> seen=new ObjectSet<String>(String.class),runs=new ObjectSet<String>(String.class),assemblyModels=new ObjectSet<String>(String.class);
		final ObjectMap<String,long[]> perAssembly=new ObjectMap<String,long[]>(String.class,long[].class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long total=0;
		try{RrnaResourceIO.header(in,SUMMARY_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==13,"Expected scan summary schema");
			final String assembly=p.parseString(0),run=p.parseString(1),model=p.parseString(2),bundle=p.parseString(3);
			require(panel.byAssembly.get(assembly)!=null && SeedOffsetTable.validSha80(bundle),"Summary assembly/bundle invalid");
			final String key=assembly+'\t'+run+'\t'+bundle+'\t'+model;require(seen.add(key) && assemblyModels.add(assembly+'\t'+model),"Each assembly/model scanned exactly once");
			final long n=p.parseLong(4);require(n>=0 && n==get(observed,key),"Summary raw count differs from retained hit rows");
			total+=n;runs.add(run);increment(perAssembly,assembly,1);scans.add(line);
		}
		}finally{require(!in.close(),"Summary read failed");}
		require(total==hits.size() && seen.size()==scans.size() && runs.size()==expectedRuns,"Hit/run denominator conservation failed");
		for(EukPanelGffValidation.Assembly a:panel.assemblies){require(get(perAssembly,a.id)==expectedModels,"Missing zero-hit or populated model scan: "+a.id);}
		for(Hit h:hits){require(seen.contains(h.summaryKey),"Hit has no declared scan denominator");}
	}
	void write(String prefix){
		assert(!scans.isEmpty()):"A proposal needs the complete scan denominator";
		RrnaResourceIO.fresh(prefix,".overlay.tsv.gz",".coverage.tsv",".counts.tsv",".scan_denominators.tsv");
		final ObjectMap<String,long[]> coverage=new ObjectMap<String,long[]>(String.class,long[].class),counts=new ObjectMap<String,long[]>(String.class,long[].class);
		final ByteStreamWriter overlay=RrnaResourceIO.writer(prefix+".overlay.tsv.gz"),hist=RrnaResourceIO.writer(prefix+".coverage.tsv"),totals=RrnaResourceIO.writer(prefix+".counts.tsv"),denom=RrnaResourceIO.writer(prefix+".scan_denominators.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			overlay.println(HIT_HEADER+"\tselection_bin\tcontig_length\tcontig_edge_contact\tfull_model_span\tmissing_model_5prime\tmissing_model_3prime\tga_pass\tsame_model_overlaps\tsame_family_other_model_overlaps\tcross_family_overlaps\tproposal_status");
			for(Hit h:hits){
				b.clear().append(h.raw).tab().append(h.bin).tab();if(h.contigLength<0){b.append("NA");}else{b.append(h.contigLength);}
				b.tab().append(h.edge).tab().append(h.full?"true":"false").tab().append(h.mf-1).tab().append(h.clen-h.mt).tab().append(h.gaPass?"true":"false")
					.tab().append(h.sameModel).tab().append(h.sameFamily).tab().append(h.crossFamily).tab().append("REVIEW_ONLY");overlay.print(b.nl());
				final String group=h.bin+'\t'+h.compartment+'\t'+h.bundle+'\t'+h.model+'\t'+h.family;
				final String shape=group+'\t'+h.clen+'\t'+h.mf+'\t'+h.mt+'\t'+h.trunc+'\t'+h.inclusion+'\t'+h.gaPass;
				increment(coverage,shape,1);
				long[] c=counts.get(group);if(c==null){c=new long[11];counts.put(group,c);}c[0]++;
				if(h.valid){c[1]++;}if(h.full){c[2]++;}if(!h.trunc.equals("no")){c[3]++;}if(h.gaPass){c[4]++;}
				if(h.edge.equals("LEFT") || h.edge.equals("RIGHT") || h.edge.equals("BOTH")){c[5]++;}
				if(h.sameModel>0){c[6]++;}if(h.sameFamily>0){c[7]++;}if(h.crossFamily>0){c[8]++;}
				if(h.compartment.equals("UNKNOWN")){c[9]++;}if(h.inclusion.equals("!")){c[10]++;}
			}
			final String groupHeader="selection_bin\tactual_compartment\tbundle_sha80\tmodel\tfamily";
			hist.println(groupHeader+"\tclen\tmodel_from\tmodel_to\ttrunc\tinclusion\tga_pass\thits");
			for(String k:sorted(coverage)){hist.print(b.clear().append(k).tab().append(coverage.get(k)[0]).nl());}
			totals.println(groupHeader+"\thits\tvalid_interval\tfull_model_span\tinfernal_truncated\tga_pass\tcontig_edge_contact\thits_with_same_model_overlap\thits_with_same_family_other_model_overlap\thits_with_cross_family_overlap\tunknown_compartment\tinfernal_included");
			long conserved=0;for(String k:sorted(counts)){b.clear().append(k);for(long n:counts.get(k)){b.tab().append(n);}totals.print(b.nl());conserved+=counts.get(k)[0];}
			require(conserved==hits.size(),"Observed-hit groups must conserve every raw row");
			denom.println("selection_bin\t"+SUMMARY_HEADER);final LineParser1 p=new LineParser1('\t');
			for(byte[] line:scans){p.set(line);denom.print(b.clear().append(bins.get(p.parseString(0))).tab().append(line).nl());}
		}finally{RrnaResourceIO.finish(overlay,hist,totals,denom);}
		System.err.println("CM_PROPOSAL_HITS_COMPLETE assemblies="+panel.assemblies.size()+" hits="+hits.size()+" scan_model_rows="+scans.size()+" physical_pairs="+pairCount+" all_hits_review_only=true no_truth_admission=true");
	}
	public static final class Hit{
		Hit(byte[] line,LineParser1 p){
			raw=line;assembly=p.parseString(0);run=p.parseString(1);lineNumber=p.parseLong(4);key=key(run,lineNumber);expectedLinks=p.parseLong(28);
			bundle=p.parseString(5);model=p.parseString(6);family=p.parseString(7);nativeSeqid=p.parseString(10);compartment=p.parseString(11);
			mf=p.parseInt(14);mt=p.parseInt(15);clen=p.parseInt(16);final long from=p.parseLong(19),to=p.parseLong(20);start=Math.min(from,to);stop=Math.max(from,to);
			strand=p.parseString(21);trunc=p.parseString(22);inclusion=p.parseString(23);final double bits=p.parseDouble(24),ga=p.parseDouble(26);
			require(lineNumber>0 && SeedOffsetTable.validSha80(bundle) && mf>0 && mt>=mf && mt<=clen && start>0,"Invalid retained hit identity/model/coordinates");
			require((strand.equals("+") && from<=to) || (strand.equals("-") && from>=to),"Retained hit strand/order mismatch");
			require(trunc.equals("no") || trunc.equals("5'") || trunc.equals("3'") || trunc.equals("5'&3'"),"Unknown truncation label");
			require(inclusion.equals("!") || inclusion.equals("?"),"Unknown inclusion marker");
			require(Double.isFinite(bits) && Double.isFinite(ga) && p.parseInt(17)==mt-mf+1,"Invalid score or model span");
			full=mf==1 && mt==clen;gaPass=bits>=ga;valid=p.termEquals("VALID_INTERVAL",13);
			summaryKey=assembly+'\t'+run+'\t'+bundle+'\t'+model;
		}
		final byte[] raw;final String assembly,run,key,bundle,model,family,nativeSeqid,compartment,strand,trunc,inclusion,summaryKey;
		final long lineNumber,start,stop,expectedLinks;final int mf,mt,clen;public final boolean full,gaPass,valid;
		String bin,edge;long contigLength=-1,sameModel,sameFamily,crossFamily,sourceLinks;
	}
	static String key(String run,long line){return run+'\t'+line;}
	static long get(ObjectMap<String,long[]> m,String k){final long[] x=m.get(k);return x==null?0:x[0];}
	static void increment(ObjectMap<String,long[]> m,String k,long n){long[] x=m.get(k);if(x==null){x=new long[1];m.put(k,x);}x[0]+=n;}
	static ArrayList<String> sorted(ObjectMap<String,long[]> m){final ArrayList<String> k=new ArrayList<String>();for(String s:m.keys()){if(s!=null){k.add(s);}}Collections.sort(k);return k;}
	public static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final String HIT_HEADER="assembly\trun_id\ttable\ttable_sha80\ttable_line\tcm_sha80\tmodel\tfamily\tsubtype\ttarget\tnative_seqid\tcompartment\tmapping_status\trange_status\tmodel_from\tmodel_to\tclen\tmodel_span\tmodel_coverage\tseq_from\tseq_to\tstrand\ttrunc\tinclusion\tbits\tevalue\tga\tga_margin\thistorical_overlaps";
	static final String PAIR_HEADER="assembly\tnative_seqid\tstrand\trun_a\tline_a\tmodel_a\trun_b\tline_b\tmodel_b\toverlap_bp\trelation";
	public static final String SUMMARY_HEADER="assembly\trun_id\tmodel\tcm_sha80\traw_hits\tvalid_interval_hits\tunresolved_hits\tout_of_range_hits\tincluded_hits\tfull_model_span_hits\tuntruncated_hits\tga_or_higher_hits\thits_with_historical_overlap";
	static final String SOURCE_HEADER=EukPanelCmEvidence.FEATURE_HEADER+"\tassociation_status\toverlapping_hits\tdistinct_models";
	static final String ASSOCIATION_HEADER="assembly\tsource_gff\tsource_line\tstatus\trun_id\ttable_line\tmodel\toverlap_bp\tsource_fraction\thit_fraction";
	EukPanelGffValidation panel=new EukPanelGffValidation();
	final ObjectMap<String,String> bins=new ObjectMap<String,String>(String.class,String.class);
	final ObjectMap<String,Hit> byHit=new ObjectMap<String,Hit>(String.class,Hit.class);
	final ObjectMap<String,long[]> observed=new ObjectMap<String,long[]>(String.class,long[].class);
	final ArrayList<Hit> hits=new ArrayList<Hit>();final ArrayList<byte[]> scans=new ArrayList<byte[]>();long pairCount;
}
