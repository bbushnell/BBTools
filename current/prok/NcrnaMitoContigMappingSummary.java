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
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.SamLine;
import structures.ByteBuilder;
import structures.LongList;

/** Query-contig evidence from ordered, single-ended, nonoverlapping tiled SAM.
 * Uses the native SAM parser. Never equates mapped tiles with whole-contig coverage
 * or reference membership with biological origin. The producer binds input hashes.
 * @author Brian, Ganyu */
public final class NcrnaMitoContigMappingSummary {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoContigMappingSummary.class,false);
		try{
			String in=null,contigs=null,out=null,names=null,controls=null,contexts=null,calls=null;int tile=500,cap=5;
			for(String arg:pp.args){final String[] p=Parse.splitOnFirst(arg,'=');require(p.length==2,"Use flag=value");
				final String key=p[0].toLowerCase(Locale.ROOT),v=p[1];
				if(key.equals("in")){in=v;}else if(key.equals("contigs")){contigs=v;}else if(key.equals("outprefix")){out=v;}
				else if(key.equals("names")){names=v;}else if(key.equals("tilelen")){tile=Parse.parseIntKMG(v);}
				else if(key.equals("controls")){controls=v;}
				else if(key.equals("contexts")){contexts=v;}else if(key.equals("calls")){calls=v;}
				else if(key.equals("maxsites")){cap=Integer.parseInt(v);}else{throw new IllegalArgumentException("Unknown argument: "+key);}}
			require(in!=null && out!=null && ((contigs!=null && contexts==null && calls==null) || (contigs==null && contexts!=null && calls!=null && names==null && controls==null)),"Bind either original contig metadata or paired context/call tables, plus SAM and fresh outputs");
			new NcrnaMitoContigMappingSummary().run(in,contigs,out,names,controls,tile,cap,contexts,calls);
		}finally{Shared.closeStream(pp.outstream);}
	}

	public void run(String sam,String metadata,String prefix,String names,int tile,int cap){
		run(sam,metadata,prefix,names,null,tile,cap);
	}
	public void run(String sam,String metadata,String prefix,String names,String controls,int tile,int cap){
		run(sam,metadata,prefix,names,controls,tile,cap,null,null);
	}
	public void run(String sam,String metadata,String prefix,String names,String controls,int tile,int cap,String contexts,String calls){
		require(contigs.isEmpty() && tile>0 && cap>0,"One run per instance with positive tile length and reported-hit cap");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);tileLength=tile;maxSites=cap;
		contextMode=contexts!=null;if(contextMode){require(calls!=null && metadata==null && names==null && controls==null,"Context mode uses only its bound context/call metadata");loadContexts(contexts,calls);}else{loadContigs(metadata,names,controls);}
		RrnaResourceIO.fresh(prefix,".alignments.tsv.gz",".tiles.tsv.gz",".contigs.tsv",".references.tsv",".summary.tsv",".calls.tsv",".call_references.tsv",".call_alignments.tsv.gz");
		final ByteStreamWriter align=RrnaResourceIO.writer(prefix+".alignments.tsv.gz"),tiles=RrnaResourceIO.writer(prefix+".tiles.tsv.gz");
		if(contextMode){callAlignments=RrnaResourceIO.writer(prefix+".call_alignments.tsv.gz");callAlignments.println("case_id\tcontext_id\ttile\treference\tflag\tprimary\tcall_matches\tcall_mismatches\tcall_insertions\tcall_internal_deletions\tcall_softclips\tcall_unknown_M\tcall_CIGAR_identity_lower\tcall_CIGAR_identity_upper\toriginal_supported_intervals0\toriginal_unknown_M_intervals0\treference_supported_intervals1\ttile_reported_YI_percent\tpasses_alignment_identity_076");}
		final ByteFile input=RrnaResourceIO.open(sam);final LineParser1 parser=new LineParser1('\t');
		final ArrayList<SamLine> group=new ArrayList<SamLine>();String previous=null;
		try{
			align.println("tile\texported_contig\ttile_start0\ttile_stop0\tflag\treference\tref_start1\tmapq\tcigar\tmatch_bases\tmismatch_bases\tinserted_bases\tdeleted_bases\tsoft_clipped_bases\tsource_supported_intervals0\treported_YI_percent\treported_raw_score\tprimary\tunknown_M_bases\tsource_unknown_M_intervals0\tCIGAR_identity_lower\tCIGAR_identity_upper"+(contextMode?"\tpasses_alignment_identity_076":""));
			tiles.println("tile\texported_contig\tstart0\tstop0\tlength\tmapped_alignments\tprimary_supported_bp\tany_reported_supported_bp\treported_references\tambiguous\treported_hit_cap_reached\tstatus\tprimary_aligned_including_unknown_M_bp\tany_aligned_including_unknown_M_bp");
			for(byte[] line;(line=input.nextLine())!=null;){if(line.length==0 || line[0]=='@'){continue;}
				parser.set(line);require(parser.terms()>=11,"SAM requires eleven mandatory fields");
				final SamLine sl=new SamLine(parser);
				if(previous!=null && !previous.equals(sl.qname)){process(group,align,tiles);group.clear();}
				previous=sl.qname;group.add(sl);
			}
			if(!group.isEmpty()){process(group,align,tiles);}
		}finally{final boolean error=input.close();if(callAlignments!=null){RrnaResourceIO.finish(callAlignments);}RrnaResourceIO.finish(align,tiles);require(!error,"SAM input read failed");}
		for(Contig c:contigs.values()){if(c!=null){require(c.observed==c.seen.length,"Missing tile(s), including unmapped/tail records: "+c.id+" observed="+c.observed+" expected="+c.seen.length);}}
		write(prefix);if(contextMode){writeCalls(prefix);}
	}
	private void loadContexts(String path,String calls){
		ByteFile input=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(input,NcrnaMitoCallContexts.CONTEXT_HEADER,path);
			for(byte[] line;(line=input.nextLine())!=null;){p.set(line);require(p.terms()==16,"Complete context metadata required");
				final int start=p.parseInt(4),stop=p.parseInt(5),length=p.parseInt(6),sourceLength=p.parseInt(8);final long acgt=p.parseLong(7),sourceAcgt=p.parseLong(9);
				require(start>=0 && stop>=start && stop<sourceLength && (long)stop-start+1==length && acgt>=0 && acgt<=length && sourceAcgt>=acgt && sourceAcgt<=sourceLength,"Context must retain valid independent source/context denominators");
				final Contig c=new Contig(p.parseString(0),p.parseString(2),p.parseString(3),length,acgt,tileLength,false);c.sourceOffset=start;c.sourceLength=sourceLength;c.sourceAcgt=sourceAcgt;c.sourceId=p.parseString(1);c.scope=p.parseString(10);
				require((c.scope.equals("WHOLE_CONTIG") && start==0 && length==sourceLength) || c.scope.equals("CALL_NEIGHBORHOOD"),"Whole-contig scope requires the actual complete source");require(contigs.put(c.id,c)==null,"Context IDs must be unique");}
		}finally{require(!input.close(),"Context metadata read failed");}
		input=RrnaResourceIO.open(calls);final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(input,NcrnaMitoCallContexts.CALL_HEADER,calls);
			for(byte[] line;(line=input.nextLine())!=null;){p.set(line);require(p.terms()==12,"Complete call/context membership required");final Contig c=contigs.get(p.parseString(1));require(c!=null,"Call refers to an unknown context");
				final Call call=new Call(p.parseString(0),c,p.parseInt(8),p.parseInt(9),p.parseString(7),p.parseString(10),p.parseString(11));
				require(seen.add(call.id) && p.termEquals(c.sourceId,2) && p.termEquals(c.assembly,3) && p.termEquals(c.original,4)
					&& p.parseInt(5)==c.sourceOffset+call.start+1 && p.parseInt(6)==c.sourceOffset+call.stop+1,"Call coordinates must bind exactly once to the original source and measured context");c.calls.add(call);}
		}finally{require(!input.close(),"Call memberships read failed");}
		require(!contigs.isEmpty(),"Context mapping cannot be empty");for(Contig c:contigs.values()){if(c!=null){require(!c.calls.isEmpty(),"Every measured context must retain its calls");}}
	}

	private void loadContigs(String path,String names,String controls){
		final ObjectSet<String> wanted=new ObjectSet<String>(String.class);
		if(names!=null){for(String name:names.split(",")){require(!name.isEmpty() && wanted.add(name),"Explicit pilot names must be unique");}}
		final ObjectSet<String> controlNames=new ObjectSet<String>(String.class);
		if(controls!=null){for(String name:controls.split(",")){require(!name.isEmpty() && controlNames.add(name) && !wanted.contains(name),"Control accession must be explicit, unique and distinct from exported IDs");}}
		int selectedControls=0;
		final ByteFile input=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{p.set(input.nextLine());require(p.terms()>=10 && p.termEquals("length",2) && p.termEquals("exported_id",6),"Expected native whole-contig evidence table");
			for(byte[] line;(line=input.nextLine())!=null;){p.set(line);require(p.terms()>=10,"Incomplete contig metadata");
				final boolean control=p.termEquals("NA",6);
				final String id=p.parseString(control?1:6);
				if(control){if(!controlNames.contains(id)){continue;}selectedControls++;}
				else if(names!=null && !wanted.contains(id)){continue;}
				final int length=p.parseInt(2);final long acgt=p.parseLong(3);
				require(length>0 && acgt>=0 && acgt<=length,"Original whole-contig denominators must be valid: "+id);
				final Contig c=new Contig(id,p.parseString(0),p.parseString(1),length,acgt,tileLength,control);
				require(contigs.put(id,c)==null,"Duplicate exported contig ID: "+id);
			}
		}finally{require(!input.close(),"Contig metadata read failed");}
		require(!contigs.isEmpty() && selectedControls==controlNames.size() && (names==null || contigs.size()==wanted.size()+selectedControls),"Every explicit pilot/control contig must exist in the bound metadata");
	}

	private void process(ArrayList<SamLine> rows,ByteStreamWriter align,ByteStreamWriter tiles){
		assert(!rows.isEmpty()):"A complete QNAME group is the attribution unit";
		final String name=rows.get(0).qname;final int underscore=name.lastIndexOf('_'),dash=name.indexOf('-',underscore+1);
		require(underscore>0 && dash>underscore+1,"Shred coordinates must survive header trimming: "+name);
		final String id=name.substring(0,underscore);final Contig c=contigs.get(id);
		require(c!=null,"Unexpected contig in SAM: "+id);
		final int start=Integer.parseInt(name.substring(underscore+1,dash)),stop=Integer.parseInt(name.substring(dash+1));
		require(start>=0 && start<c.length && start%tileLength==0 && stop==(int)Math.min((long)start+tileLength,c.length)-1,
			"Require the declared complete nonoverlapping tiling, including its short tail: "+name);
		final int index=start/tileLength,len=stop-start+1;
		require(!c.seen[index],"Repeated/noncontiguous SAM tile group: "+name);c.seen[index]=true;c.observed++;
		allIntervals.clear();primaryIntervals.clear();allPossible.clear();primaryPossible.clear();rejectedIntervals.clear();touched.clear();int primary=0,mapped=0;boolean ambiguous=false;
		for(SamLine sl:rows){
			require(sl.qname.equals(name) && (sl.flag&1)==0 && !sl.supplementary(),"Only grouped single-end mappings without supplementary clipping are supported");
			require(sl.seq==null || sl.seq.length==len,"SAM sequence length must agree with the original tile coordinates: "+name);
			final boolean isPrimary=sl.primary();if(isPrimary){primary++;}
			if(!sl.mapped()){
				require(isPrimary && rows.size()==1,"Unmapped tile must have exactly one primary SAM record");
				continue;
			}
			require(sl.pos>0 && sl.hasCigar() && sl.rnameS()!=null,"Mapped SAM row requires reference, position and CIGAR");
			mapped++;ambiguous|=sl.ambiguous();
			parseCigar(sl.cigar,start,len,sl.strand()==1);
			final boolean passes=!contextMode || passesIdentity(matches,mismatches,insertions,deletions,unknown);
			if(passes){allIntervals.addAll(intervals);allPossible.addAll(intervals);allPossible.addAll(unknownIntervals);
				if(isPrimary){primaryIntervals.addAll(intervals);primaryPossible.addAll(intervals);primaryPossible.addAll(unknownIntervals);c.unknown+=unknown;c.matches+=matches;c.mismatches+=mismatches;c.insertions+=insertions;c.deletions+=deletions;c.clipped+=clipped;}
			}else{rejectedIntervals.addAll(intervals);}
			final String ref=sl.rnameS();Reference support=c.references.get(ref);
			if(support==null){support=new Reference(ref);c.references.put(ref,support);}
			if(support.epoch!=c.observed){support.epoch=c.observed;support.intervals.clear();support.possible.clear();support.rejected.clear();touched.add(support);}
			if(passes){support.intervals.addAll(intervals);support.possible.addAll(intervals);support.possible.addAll(unknownIntervals);support.passingAlignments++;}else{support.rejected.addAll(intervals);}
			support.alignments++;if(sl.strand()==1){support.reverse++;}else{support.forward++;}
			b.clear().append(name).tab().append(id).tab().append(start).tab().append(stop).tab().append(sl.flag)
				.tab().append(ref).tab().append(sl.pos).tab().append(sl.mapq).tab().append(sl.cigar)
				.tab().append(matches).tab().append(mismatches).tab().append(insertions).tab().append(deletions).tab().append(clipped).tab();
			appendIntervals(b,intervals);b.tab().append(tag(sl,"YI:f:")).tab().append(tag(sl,"YR:i:")).tab().append(isPrimary).tab().append(unknown).tab();
			appendIntervals(b,unknownIntervals);appendBounds(b,matches,mismatches,insertions,deletions,unknown);if(contextMode){b.tab().append(passes);}align.print(b.nl());
			for(Call call:c.calls){if(call.start<=stop && call.stop>=start){observeCall(call,sl,start,len,isPrimary,passes);}}
		}
		require(primary==1,"Exactly one primary per input tile is required to avoid identity double counting: "+name+" primaries="+primary);
		final long pc=union(primaryIntervals),any=union(allIntervals);require(pc<=any && any<=len,"Reported alignments cannot cover more source bases than one tile");
		final long pp=union(primaryPossible),ap=union(allPossible);require(pc<=pp && any<=ap && pp<=ap && ap<=len,"Unknown M support is only a possible aligned-base upper bound");
		for(Reference ref:touched){ref.covered+=union(ref.intervals);ref.possibleCovered+=union(ref.possible);ref.rejectedCovered+=union(ref.rejected);ref.tiles++;}
		c.rejectedCovered+=union(rejectedIntervals);
		final boolean cap=mapped>=maxSites;ambiguous|=mapped>1;
		c.primaryCovered+=pc;c.anyCovered+=any;c.primaryPossible+=pp;c.anyPossible+=ap;c.mappedTiles+=mapped>0?1:0;c.ambiguousTiles+=ambiguous?1:0;c.cappedTiles+=cap?1:0;c.alignments+=mapped;
		tiles.print(b.clear().append(name).tab().append(id).tab().append(start).tab().append(stop).tab().append(len)
			.tab().append(mapped).tab().append(pc).tab().append(any).tab().append(touched.size()).tab().append(ambiguous)
			.tab().append(cap).tab().append(mapped==0?"UNMAPPED":"MAPPED_REPORTED_HITS_ONLY").tab().append(pp).tab().append(ap).nl());
	}
	/** Restrict primary edit counts to a call, not the whole overlapping tile.
	 * A deletion belongs to the call only when BOTH flanking query bases are inside it. */
	private void observeCall(Call call,SamLine sl,int tileStart,int length,boolean primary,boolean passes){
		final boolean reverse=sl.strand()==1;final int left=Math.max(call.start,tileStart)-tileStart,right=Math.min(call.stop,tileStart+length-1)-tileStart;
		final int qlo=reverse?length-1-right:left,qhi=reverse?length-1-left:right;
		require(qlo>=0 && qhi>=qlo && qhi<length,"Call/tile intersection must be a nonempty query interval");
		callKnown.clear();callUnknown.clear();callRef.clear();long mat=0,sub=0,ins=0,del=0,clip=0,unk=0;int q=0,n=0;long ref=sl.pos;
		for(int i=0;i<sl.cigar.length();i++){
			final char op=sl.cigar.charAt(i);if(op>='0' && op<='9'){n=n*10+op-'0';continue;}
			if(op=='D'){if(q>qlo && q<=qhi){del+=n;}ref+=n;}
			else{final int lo=Math.max(q,qlo),hi=Math.min(q+n-1,qhi),count=Math.max(0,hi-lo+1);
				if(op=='=' || op=='X' || op=='M'){
					if(count>0){final int sourceStart=call.contig.sourceOffset+tileStart+(reverse?length-1-hi:lo),sourceStop=sourceStart+count-1;
						(op=='M'?callUnknown:callKnown).add(pack(sourceStart,sourceStop));
						if(op!='M'){require(ref+hi-q<=Integer.MAX_VALUE,"Reference coordinate cannot exceed packed interval range");callRef.add(pack((int)(ref+lo-q),(int)(ref+hi-q)));}
						if(op=='='){mat+=count;}else if(op=='X'){sub+=count;}else{unk+=count;}}
					ref+=n;
				}else if(op=='I'){ins+=count;}else if(op=='S'){clip+=count;}
				else{throw new IllegalArgumentException("Unsupported call-level CIGAR: "+sl.cigar);}q+=n;
			}n=0;
		}
		require(n==0 && q==length,"Previously validated CIGAR must conserve query during call restriction");
		if(passes){call.known.addAll(callKnown);call.possible.addAll(callKnown);call.possible.addAll(callUnknown);
			if(primary){call.primaryKnown.addAll(callKnown);call.matches+=mat;call.mismatches+=sub;call.insertions+=ins;call.deletions+=del;call.clipped+=clip;call.unknown+=unk;}
		}else{call.rejectedKnown.addAll(callKnown);}
		final String name=sl.rnameS();CallReference r=call.references.get(name);if(r==null){r=new CallReference();call.references.put(name,r);}r.alignments++;
		if(passes){r.known.addAll(callKnown);r.refIntervals.addAll(callRef);r.passingAlignments++;if(primary){r.primaryAlignments++;}}else{r.rejectedKnown.addAll(callKnown);}
		b.clear().append(call.id).tab().append(call.contig.id).tab().append(sl.qname).tab().append(name).tab().append(sl.flag).tab().append(primary)
			.tab().append(mat).tab().append(sub).tab().append(ins).tab().append(del).tab().append(clip).tab().append(unk);appendBounds(b,mat,sub,ins,del,unk);b.tab();appendIntervals(b,callKnown);b.tab();appendIntervals(b,callUnknown);b.tab();appendIntervals(b,callRef);
		callAlignments.print(b.tab().append(tag(sl,"YI:f:")).tab().append(passes).nl());
	}
	private static long pack(int start,int stop){require(start>=0 && stop>=start,"Nonnegative inclusive source/reference intervals required");return ((long)start<<32)|(stop&0xffffffffL);}
	private void writeCalls(String prefix){
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".calls.tsv"),refs=RrnaResourceIO.writer(prefix+".call_references.tsv");int count=0;
		try{out.println("case_id\tcontext_id\texported_contig\tassembly\toriginal_contig\tcall_start1\tcall_stop1\tstrand\tcall_length\tprimary_supported_bp\tany_supported_bp\tany_possible_including_M_bp\tcall_coverage\tprimary_matches\tprimary_mismatches\tprimary_insertions\tprimary_internal_deletions\tprimary_unknown_M\tprimary_identity\tidentity_lower\tidentity_upper\tsource_length\tsource_ACGT\tcontext_length\tcontext_ACGT\tcontext_known_coverage\tscope\torigin\taudit_id\twinning_model\tbelow_identity_supported_call_bp\talignment_identity_min");
			refs.println("case_id\tcontext_id\treference\tknown_supported_call_bp\tcall_coverage\treported_alignments\tprimary_alignments\treference_supported_intervals1\tbelow_identity_supported_call_bp\tpassing_alignments");
			for(String id:keys(contigs)){final Contig c=contigs.get(id);for(Call call:c.calls){count++;final int len=call.stop-call.start+1;final long known=union(call.known),primary=union(call.primaryKnown),possible=union(call.possible);require(primary<=known && known<=possible && possible<=len,"Call support unions cannot exceed the exact call denominator");
				final double coverage=known/(double)len,contextCoverage=c.anyCovered/(double)c.length;final String origin=classify(c.sourceLength,c.scope.equals("WHOLE_CONTIG"),contextCoverage,coverage);
				b.clear().append(call.id).tab().append(c.id).tab().append(c.sourceId).tab().append(c.assembly).tab().append(c.original).tab().append(c.sourceOffset+call.start+1).tab().append(c.sourceOffset+call.stop+1).tab().append(call.strand).tab().append(len)
					.tab().append(primary).tab().append(known).tab().append(possible).tab().append(coverage,9).tab().append(call.matches).tab().append(call.mismatches).tab().append(call.insertions).tab().append(call.deletions).tab().append(call.unknown).tab();
				final long denom=call.matches+call.mismatches+call.insertions+call.deletions;if(denom==0 || call.unknown>0){b.append("NA");}else{b.append(call.matches/(double)denom,9);}appendBounds(b,call.matches,call.mismatches,call.insertions,call.deletions,call.unknown);
				out.print(b.tab().append(c.sourceLength).tab().append(c.sourceAcgt).tab().append(c.length).tab().append(c.acgt).tab().append(contextCoverage,9).tab().append(c.scope).tab().append(origin).tab().append(call.audit).tab().append(call.model).tab().append(union(call.rejectedKnown)).append("\t0.76\n"));
				for(String ref:keys(call.references)){final CallReference r=call.references.get(ref);final long supported=union(r.known);b.clear().append(call.id).tab().append(c.id).tab().append(ref).tab().append(supported).tab().append(supported/(double)len,9).tab().append(r.alignments).tab().append(r.primaryAlignments).tab();appendIntervals(b,r.refIntervals);refs.print(b.tab().append(union(r.rejectedKnown)).tab().append(r.passingAlignments).nl());}
			}}
		}finally{RrnaResourceIO.finish(out,refs);}
		System.err.println("MITO_CALL_MAPPING_PASS calls="+count+" contexts="+contigs.size()+" coverage=known_query_union origin_descriptive_only=true");
	}
	/** G11 October22: descriptive origin bins, never a functionality determination. */
	public static String classify(int sourceLength,boolean whole,double wholeCoverage,double callCoverage){
		require(sourceLength>0 && wholeCoverage>=0 && wholeCoverage<=1 && callCoverage>=0 && callCoverage<=1,"Origin bins require valid measured coverage");
		if(callCoverage<.20){return "UNMAPPED";}if(callCoverage>=.80){if(sourceLength>=200000){return "NUMT";}require(whole,"Small-contig origin requires actual whole-contig coverage");return wholeCoverage>=.80?"MITO_CONTIG":"NUMT";}return "AMBIGUOUS";
	}
	/** Conservative whole-alignment identity; an ambiguous M cannot supply a match. */
	public static boolean passesIdentity(long matches,long mismatch,long ins,long del,long unknown){
		require(matches>=0 && mismatch>=0 && ins>=0 && del>=0 && unknown>=0,"Identity requires nonnegative actual CIGAR counts");final long denom=matches+mismatch+ins+del+unknown;return denom>0 && matches/(double)denom>=.76;
	}

	/** Convert =/X query positions to original forward-contig coordinates. */
	private void parseCigar(String cigar,int start,int length,boolean reverse){
		intervals.clear();unknownIntervals.clear();matches=mismatches=insertions=deletions=clipped=unknown=0;int q=0,n=0;
		for(int i=0;i<cigar.length();i++){
			final char op=cigar.charAt(i);if(op>='0' && op<='9'){require(n<=(Integer.MAX_VALUE-(op-'0'))/10,"CIGAR length overflow");n=n*10+(op-'0');continue;}
			require(n>0,"Each CIGAR operation must have positive length: "+cigar);
			if(op=='=' || op=='X' || op=='M'){
				require(n<=length-q,"CIGAR consumes beyond original tile");
				final int lo=start+(reverse?length-q-n:q),hi=lo+n-1;
				(op=='M'?unknownIntervals:intervals).add(((long)lo<<32)|(hi&0xffffffffL));
				if(op=='='){matches+=n;}else if(op=='X'){mismatches+=n;}else{unknown+=n;}q+=n;
			}else if(op=='I' || op=='S'){
				require(n<=length-q,"Insertion/clip consumes beyond original tile");q+=n;if(op=='I'){insertions+=n;}else{clipped+=n;}
			}else if(op=='D'){deletions+=n;}
			else{throw new IllegalArgumentException("Unsupported hard clip/splice/padding CIGAR operation: "+cigar);}
			n=0;
		}
		require(n==0 && q==length,"CIGAR must account for every original tile base: "+cigar+" query="+q+" expected="+length);
	}

	public static long union(LongList values){
		assert(values!=null):"Coverage requires explicit aligned query intervals";
		if(values.size==0){return 0;}values.sort();long total=0;int lo=(int)(values.array[0]>>>32),hi=(int)values.array[0];
		for(int i=1;i<values.size;i++){final int a=(int)(values.array[i]>>>32),z=(int)values.array[i];
			if((long)a<=1L+hi){hi=Math.max(hi,z);}else{total+=(long)hi-lo+1;lo=a;hi=z;}}
		return total+(long)hi-lo+1;
	}
	private static void appendIntervals(ByteBuilder b,LongList values){
		assert(values!=null):"Alignment rows retain the actual source intervals before aggregation";
		for(int i=0;i<values.size;i++){if(i>0){b.append(',');}b.append((int)(values.array[i]>>>32)).append('-').append((int)values.array[i]);}
	}
	private static String tag(SamLine row,String prefix){
		assert(row!=null && prefix!=null):"Read native optional tags without interpreting absent values";
		if(row.optional!=null){for(String value:row.optional){if(value.startsWith(prefix)){return value.substring(prefix.length());}}}return "NA";
	}
	private static void appendBounds(ByteBuilder b,long matches,long mismatch,long ins,long del,long unknown){
		final long denom=matches+mismatch+ins+del+unknown;assert(denom>=0):"Nonnegative CIGAR evidence counts required";
		if(denom==0){b.append("\tNA\tNA");}else{b.tab().append(matches/(double)denom,9).tab().append((matches+unknown)/(double)denom,9);}
	}
	private void write(String prefix){
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".contigs.tsv"),refs=RrnaResourceIO.writer(prefix+".references.tsv"),summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		long bases=0,tiles=0,mapped=0;
		try{
			out.println("exported_id\tassembly\toriginal_contig\tlength\tACGT\ttiles\tmapped_tiles\tprimary_supported_bp\tany_reported_supported_bp\t"+(contextMode?"primary_fraction_context\tany_fraction_context":"primary_fraction_full_contig\tany_fraction_full_contig")+"\tprimary_CIGAR_matches\tprimary_CIGAR_mismatches\tprimary_insertions\tprimary_deletions\tprimary_softclips\tprimary_CIGAR_match_fraction\tambiguous_tiles\treported_hit_cap_tiles\treported_alignments\tinput_role\tprimary_unknown_M_bases\tprimary_CIGAR_identity_lower\tprimary_CIGAR_identity_upper\tprimary_aligned_including_M_bp\tany_aligned_including_M_bp"+(contextMode?"\tsource_exported_id\tsource_start0\tsource_length\tsource_ACGT\tscope\tbelow_identity_supported_bp":""));
			refs.println("exported_id\treference\treported_supported_query_bp\t"+(contextMode?"fraction_context":"fraction_full_contig")+"\ttiles\treported_alignments\tforward_alignments\treverse_alignments\taligned_including_unknown_M_bp"+(contextMode?"\tbelow_identity_supported_bp\tpassing_alignments":""));
			for(String id:keys(contigs)){final Contig c=contigs.get(id);bases+=c.length;tiles+=c.observed;mapped+=c.mappedTiles;
				b.clear().append(id).tab().append(c.assembly).tab().append(c.original).tab().append(c.length).tab().append(c.acgt)
					.tab().append(c.observed).tab().append(c.mappedTiles).tab().append(c.primaryCovered).tab().append(c.anyCovered)
					.tab().append(c.primaryCovered/(double)c.length,9).tab().append(c.anyCovered/(double)c.length,9)
					.tab().append(c.matches).tab().append(c.mismatches).tab().append(c.insertions).tab().append(c.deletions).tab().append(c.clipped).tab();
				final long denom=c.matches+c.mismatches+c.insertions+c.deletions;if(denom==0 || c.unknown>0){b.append("NA");}else{b.append(c.matches/(double)denom,9);}
				b.tab().append(c.ambiguousTiles).tab().append(c.cappedTiles).tab().append(c.alignments).tab().append(c.control?"EXPLICIT_REFERENCE_CONTROL":"UNRESOLVED_CANDIDATE_CONTIG").tab().append(c.unknown);
				appendBounds(b,c.matches,c.mismatches,c.insertions,c.deletions,c.unknown);b.tab().append(c.primaryPossible).tab().append(c.anyPossible);
				if(contextMode){b.tab().append(c.sourceId).tab().append(c.sourceOffset).tab().append(c.sourceLength).tab().append(c.sourceAcgt).tab().append(c.scope).tab().append(c.rejectedCovered);}out.print(b.nl());
				for(String name:keys(c.references)){final Reference r=c.references.get(name);b.clear().append(id).tab().append(name).tab().append(r.covered)
					.tab().append(r.covered/(double)c.length,9).tab().append(r.tiles).tab().append(r.alignments).tab().append(r.forward).tab().append(r.reverse).tab().append(r.possibleCovered);
					if(contextMode){b.tab().append(r.rejectedCovered).tab().append(r.passingAlignments);}refs.print(b.nl());}
			}
			summary.println("contigs\tquery_bases\ttiles\tmapped_tiles\ttile_length\treported_hit_cap\torigin_classification\tinterpretation");
			summary.print(b.clear().append(contigs.size()).tab().append(bases).tab().append(tiles).tab().append(mapped).tab().append(tileLength).tab().append(maxSites)
				.append("\tUNRESOLVED\tReported alignment coverage; not exhaustive homology, ANI or a single collinear whole-contig alignment\n"));
		}finally{RrnaResourceIO.finish(out,refs,summary);}
		System.err.println("MITO_CONTIG_MAPPING_SUMMARY_PASS contigs="+contigs.size()+" bases="+bases+" tiles="+tiles+" mapped="+mapped);
	}
	private static <V> ArrayList<String> keys(ObjectMap<String,V> map){final ArrayList<String> keys=new ArrayList<String>();for(String k:map.keys()){if(k!=null){keys.add(k);}}Collections.sort(keys);return keys;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Contig{
		Contig(String id_,String a,String o,int len,long acgt_,int tile,boolean control_){id=id_;assembly=a;original=o;length=len;acgt=acgt_;control=control_;seen=new boolean[(int)(((long)len+tile-1)/tile)];sourceId=id;sourceLength=len;sourceAcgt=acgt;}
		final String id,assembly,original;final int length;final long acgt;final boolean control;final boolean[] seen;int observed;
		String sourceId,scope="WHOLE_CONTIG";int sourceOffset,sourceLength;long sourceAcgt;final ArrayList<Call> calls=new ArrayList<Call>();
		long primaryCovered,anyCovered,primaryPossible,anyPossible,unknown,mappedTiles,ambiguousTiles,cappedTiles,alignments,matches,mismatches,insertions,deletions,clipped,rejectedCovered;
		final ObjectMap<String,Reference> references=new ObjectMap<String,Reference>(String.class,Reference.class);
	}
	private static final class Reference{
		Reference(String name){id=name;}final String id;int epoch;long covered,possibleCovered,tiles,alignments,forward,reverse,rejectedCovered,passingAlignments;final LongList intervals=new LongList(),possible=new LongList(),rejected=new LongList();
	}
	private static final class Call{
		Call(String i,Contig c,int a,int z,String s,String audit_,String model_){id=i;contig=c;start=a;stop=z;strand=s;audit=audit_;model=model_;require(a>=0 && z>=a && z<c.length && (s.equals("+") || s.equals("-")),"Each call must fit its context with an explicit strand");}
		final String id,strand,audit,model;final Contig contig;final int start,stop;long matches,mismatches,insertions,deletions,clipped,unknown;
		final LongList known=new LongList(),primaryKnown=new LongList(),possible=new LongList(),rejectedKnown=new LongList();final ObjectMap<String,CallReference> references=new ObjectMap<String,CallReference>(String.class,CallReference.class);
	}
	private static final class CallReference{long alignments,primaryAlignments,passingAlignments;final LongList known=new LongList(),refIntervals=new LongList(),rejectedKnown=new LongList();}
	private boolean contextMode;private ByteStreamWriter callAlignments;private final LongList callKnown=new LongList(),callUnknown=new LongList(),callRef=new LongList();
	private int tileLength,maxSites;private long matches,mismatches,insertions,deletions,clipped,unknown;
	private final ObjectMap<String,Contig> contigs=new ObjectMap<String,Contig>(String.class,Contig.class);
	private final LongList intervals=new LongList(),unknownIntervals=new LongList(),allIntervals=new LongList(),primaryIntervals=new LongList(),allPossible=new LongList(),primaryPossible=new LongList(),rejectedIntervals=new LongList();
	private final ArrayList<Reference> touched=new ArrayList<Reference>();private final ByteBuilder b=new ByteBuilder();
}
