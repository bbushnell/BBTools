package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import structures.ByteBuilder;
import shared.Shared;

/** Reuses dispositions only at identical original loci, projects explicitly bound
 * independent hits onto added calls, and reports production precision without
 * adding TPs. Historical GA support is not presented as three new model scans.
 * @author Brian Bushnell, Ganyu
 */
public final class CobalaminProductionPrecision {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminProductionPrecision.class, false);
		final ObjectMap<String,String> a=new ObjectMap<String,String>(String.class, String.class);
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT);
			require((" old new review sources decisions calls summary context features neighborhood outprefix oldcount newcount fpcount ").contains(" "+key+" ") && !a.contains(key), "Unknown/repeated argument: "+key);a.put(key, arg.substring(eq+1));}
		for(String key:new String[]{"old","new","review","sources","decisions","calls","summary","context","features","neighborhood","outprefix","oldcount","newcount","fpcount"}){require(a.get(key)!=null, "Missing argument: "+key);}
		Shared.setThreads(1);
		try{new CobalaminProductionPrecision().run(a);}finally{Shared.closeStream(pp.outstream);}
	}
	private void run(ObjectMap<String,String> a){
		final ArrayList<Locus> old=metadata(a.get("old")), next=metadata(a.get("new"));loci=next;
		require(old.size()==Integer.parseInt(a.get("oldcount")) && next.size()==Integer.parseInt(a.get("newcount")) && next.size()>old.size(), "Bind the original and augmented exact-locus populations");
		final ObjectMap<String,Locus> oldIds=new ObjectMap<String,Locus>(String.class, Locus.class);
		for(Locus l:next){require(!byId.contains(l.id) && !byKey.contains(l.key()), "Duplicate augmented locus identity");byId.put(l.id,l);byKey.put(l.key(),l);}
		for(Locus p:old){
			final Locus l=byKey.get(p.key());require(l!=null && sameMetadata(p,l) && l.oldId==null && !oldIds.contains(p.id), "Every old locus must map exactly once with unchanged metadata, never by row number");
			l.oldId=p.id;oldIds.put(p.id,l);
		}
		read(a.get("review"), REVIEW_HEADER, p->{final Locus l=oldIds.get(p.parseString(0));
			require(l!=null && l.category==null && category(p.parseString(1)), "Old dispositions must cover each exact original locus once");
			l.category=p.parseString(1);l.reason=p.parseString(2);
		});
		for(Locus l:next){if(l.oldId!=null){require(l.category!=null, "Missing original disposition");}}
		context(a.get("context"));features(a.get("features"));neighborhood(a.get("neighborhood"));
		read(a.get("decisions"), CobalaminPrecisionReview.DECISION_HEADER, p->{final Locus l=byId.get(p.parseString(0));
			require(l!=null && l.oldId==null && l.decision==null && p.termEquals(GEOMETRY,1) && !p.termEquals("",2), "Only explicit added-locus annotated-geometry decisions are supported");
			l.decision=p.parseString(1);l.reason=p.parseString(2);
		});
		read(a.get("sources"), SOURCE_HEADER, p->source(p));
		for(Locus l:next){if(l.oldId!=null){continue;}
			require(l.sources>0, "Added locus has no independently bound scan: "+l.id);
			if(l.decision!=null){require(l.fullContained>0 && l.annotatedContained, "Annotated geometry requires a containing full-GA hit AND same-strand source annotation");l.category=GEOMETRY;}
			else{require(l.full==0 && l.scanMask==63, "New full-GA support needs a decision; unsupported calls need all three DEFAULT/MAX scans");
				l.category=l.weak>0?PARTIAL:NONE;l.reason=l.weak>0?"Overlapping below-GA/partial evidence; no optimistic exclusion":"No overlapping reported hit in six scans; not proof of biological absence";}
		}
		calls(a.get("calls"), Integer.parseInt(a.get("fpcount")));summaries(a.get("summary"));write(a.get("outprefix"));
	}
	private static ArrayList<Locus> metadata(String path){
		final ArrayList<Locus> list=new ArrayList<Locus>();final ObjectSet<String> ids=new ObjectSet<String>(String.class), keys=new ObjectSet<String>(String.class);
		read(path,NcrnaReviewUnion.LOCUS_HEADER,p->{final Locus l=new Locus(p);require(ids.add(l.id) && keys.add(l.key()),"Source metadata must have unique IDs and exact loci");list.add(l);});return list;
	}
	private void context(String path){
		read(path,NcrnaReviewJoin.CONTEXT_HEADER,p->{final Locus l=locus(p);
			require(l.context==null && p.termEquals(l.genome,1) && p.termEquals(l.assembly,2) && p.termEquals(l.seq,3) && p.parseInt(4)==l.start && p.parseInt(5)==l.stop && p.termEquals(l.strand,6) && p.parseInt(7)==l.length() && p.termEquals(l.fasta,24),"Feature context must retain exact augmented identity");
			l.context=p.line().clone();l.gff=p.parseString(22);l.gffHash=p.parseString(23);l.fastaHash=p.parseString(25);l.expectedFeatures=p.parseInt(11);l.expectedRegions=p.parseInt(12);
			require(l.gffHash.matches("[0-9a-f]{20}") && l.fastaHash.matches("[0-9a-f]{20}") && l.expectedFeatures>=0 && l.expectedRegions>=0,"Context requires source hashes and nonnegative counts");
		});for(Locus l:loci){require(l.context!=null,"Missing source context: "+l.id);}
	}
	private void features(String path){
		read(path,NcrnaReviewJoin.OVERLAP_HEADER,p->{final Locus l=locus(p);final String kind=p.parseString(8), strand=p.parseString(19);
			final int start=p.parseInt(16), stop=p.parseInt(17), bp=p.parseInt(10);
			require(p.termEquals(l.genome,1) && p.termEquals(l.assembly,2) && p.termEquals(l.fasta,3) && p.termEquals(l.fastaHash,4) && p.termEquals(l.gff,5) && p.termEquals(l.gffHash,6) && p.termEquals(l.seq,13),"Feature sources must bind the same genome");
			require(start>0 && stop>=start && p.parseInt(7)>0 && l.featureLines.add(p.parseString(7)) && bp>0 && bp==NcrnaReviewFeatureContext.overlap(l.start,l.stop,start,stop,l.contigLength),"Feature geometry or row conservation changed");
			final String relation=strand.equals("+") || strand.equals("-")?(strand.equals(l.strand)?"SAME":"OPPOSITE"):"UNKNOWN";
			require(p.termEquals(relation,11),"Feature strand relation differs");
			if(kind.equals("REGION_METADATA")){l.regions++;}else{l.features++;}
			final String attrs=p.parseString(21);
			if(!kind.equals("REGION_METADATA") && !kind.equals("PARENT_ANNOTATION") && relation.equals("SAME") && start<=l.start && stop>=l.stop && b12(attrs)){l.annotatedContained=true;}
			l.featureRows.add(p.line().clone());
		});for(Locus l:loci){require(l.features==l.expectedFeatures && l.regions==l.expectedRegions,"Feature detail must reproduce source context counters");}
	}
	private void neighborhood(String path){
		read(path,CobalaminGenomicContext.HEADER,p->{final Locus l=locus(p);require(l.neighborhood==null,"Repeated neighborhood would multiply calls");l.neighborhood=p.line().clone();});
		for(Locus l:loci){require(l.neighborhood!=null,"Missing native neighborhood: "+l.id);}
	}
	private void source(LineParser1 p){
		final Locus l=byId.get(p.parseString(0));final String sourceId=p.parseString(1), format=p.parseString(2), metadata=p.parseString(3), hits=p.parseString(4), model=p.parseString(5), stage=p.parseString(6), encoding=p.parseString(7);
		final int m=model(model), s=stage.equals("DEFAULT")?0:stage.equals("MAX")?1:stage.equals("HISTORICAL_GA")?2:-1;
		require(l!=null && l.oldId==null && m>=0 && s>=0 && (encoding.equals("RAW") || encoding.equals("PROJECTED")) && sourceKeys.add(l.id+'\t'+model+'\t'+stage),"Added scan bindings must be explicit and unique per model/stage");
		final Locus[] selected={null};
		if(format.equals("REVIEW")){for(Locus c:metadata(metadata)){if(c.id.equals(sourceId)){selected[0]=c;}}}
		else{require(format.equals("PILOT") && encoding.equals("RAW") && s==2 && m==0,"Historical pilot means the explicit RF00174 GA scan");
			read(metadata,NcrnaCmLocusAudit.HEADER,q->{if(q.termEquals(sourceId,0)){require(selected[0]==null,"Duplicate historical source locus");selected[0]=Locus.pilot(q);}});
		}
		final Locus c=selected[0];
		require(c!=null && c.genome.equals(l.genome) && c.assembly.equals(l.assembly) && c.seq.equals(l.seq) && c.strand.equals(l.strand) && c.contigLength==l.contigLength && l.start>=c.contextStart() && l.stop<=c.contextStop(),"Reuse requires actual same-genome/strand scanned-window containment, not proximity");
		final int initial=projected.size();
		if(encoding.equals("RAW")){
			final RrnaCmTruthBuilder.Criteria criteria=new RrnaCmTruthBuilder.Criteria(model,"cobalamin",CLEN[m],CLEN[m],GA[m]);
			for(RrnaCmTruthBuilder.Hit h:RrnaCmTruthBuilder.readHits(new RrnaCmTruthBuilder.Genome("review","contexts",hits,null),criteria)){
				if(!h.seqid.equals(sourceId)){continue;}
				final int lo=Math.min(h.seqFrom,h.seqTo), hi=Math.max(h.seqFrom,h.seqTo);
				require(lo>0 && hi<=c.contextLength,"Raw CM hit lies outside its actual scanned sequence");
				final int start=c.strand.equals("+")?c.contextStart()+lo-1:c.contextStop()-hi+1;
				final int stop=c.strand.equals("+")?c.contextStart()+hi-1:c.contextStop()-lo+1;
				final String strand=h.strand=='+'?c.strand:c.strand.equals("+")?"-":"+";
				addHit(l,c,hits,h.line,m,stage,start,stop,strand,h.mdlFrom,h.mdlTo,h.trunc,h.score,h.evalue,h.included);
			}
		}else{
			require(format.equals("REVIEW") && s<2,"Projected review evidence must retain its declared DEFAULT/MAX stage");
			read(hits,NcrnaReviewJoin.HIT_HEADER,q->{if(!q.termEquals(sourceId,0)){return;}
				require(q.termEquals(c.seq,2),"Projected CM contig mismatch");
				final int start=q.parseInt(3), stop=q.parseInt(4);final String strand=q.parseString(5);
				final int oldOverlap=strand.equals(c.strand)?overlap(start,stop,c.start,c.stop):0;
				require(q.parseInt(12)==oldOverlap,"Stored raw-hit projection must reproduce its original overlap");
				addHit(l,c,hits,q.parseInt(1),m,stage,start,stop,strand,q.parseInt(6),q.parseInt(7),q.parseString(8),q.parseString(9),q.parseString(10),bool(q,11));
			});
		}
		l.sources++;if(s<2){l.scanMask|=1<<(s*3+m);}
		sourceRows.add(new ByteBuilder().append(l.id).tab().append(sourceId).tab().append(format).tab().append(metadata).tab().append(hits).tab().append(model).tab().append(stage).tab().append(c.contextStart()).tab().append(c.contextStop()).tab().append(projected.size()-initial).toBytes());
	}
	private void addHit(Locus l,Locus c,String file,int line,int m,String stage,int start,int stop,String strand,int from,int to,String trunc,String bitsText,String evalue,boolean included){
		final double bits=Double.parseDouble(bitsText), e=Double.parseDouble(evalue);
		require(line>0 && start>=c.contextStart() && stop>=start && stop<=c.contextStop() && (strand.equals("+") || strand.equals("-")) && from>=1 && to>=from && to<=CLEN[m] && Double.isFinite(bits) && Double.isFinite(e) && e>=0,"Projected hit must retain valid source geometry/model/score");
		final int bp=strand.equals(l.strand)?overlap(start,stop,l.start,l.stop):0;
		final boolean full=bp>0 && included && bits>=GA[m] && from==1 && to==CLEN[m] && trunc.equals("no");
		if(full){l.full++;if(start<=l.start && stop>=l.stop){l.fullContained++;}}
		if(bp>0 && !full){l.weak++;}
		projected.add(new ByteBuilder().append(l.id).tab().append(c.id).tab().append(file).tab().append(line).tab().append(MODELS[m]).tab().append(stage).tab().append(l.seq).tab().append(start).tab().append(stop).tab().append(strand).tab().append(from).tab().append(to).tab().append(trunc).tab().append(bitsText).tab().append(evalue).tab().append(included).tab().append(bp).tab().append(full).toBytes());
	}
	private void calls(String path,int expectedFp){
		read(path,NcrnaReviewUnion.PANEL_CALL_HEADER,p->{if(!p.termEquals("FULL_CM",0)){return;}
			final String scope=p.parseString(1), id=p.parseString(2);final int s=scope.equals("union")?3:model(scope);require(s>=0,"Unknown native grading scope");
			final Call c=new Call(p);require(!calls[s].contains(id),"Repeated native call ID within scope");calls[s].put(id,c);ordered[s].add(c);
			if(c.outcome.equals("APPARENT_FP") && s==3){final Locus l=byKey.get(c.key());require(l!=null && c.genome.equals(l.genome) && l.production==null,"Every production FP must bind exactly one augmented locus");l.production=c;c.locus=l;fps++;}
		});require(fps==expectedFp,"Production FP denominator changed");
		int partition=0;for(int s=0;s<3;s++){for(Call c:ordered[s]){final Call union=calls[3].get(c.id);require(union!=null && c.same(union) && model(c.family)==s,"Subtype calls must partition the exact union geometry/model");c.locus=union.locus;partition++;}}
		require(partition==ordered[3].size(),"Native subtype partitions must conserve every union call");
	}
	private void summaries(String path){
		read(path,SUMMARY_HEADER,p->{if(!p.termEquals("FULL_CM",0)){return;}final String scope=p.parseString(1);final int s=scope.equals("union")?3:model(scope);
			require(s>=0 && !summarySeen[s],"One native summary per scope is required");summarySeen[s]=true;
			int tp=0,fp=0,unknown=0;for(Call c:ordered[s]){if(c.outcome.equals("TP")){tp++;}else if(c.outcome.equals("APPARENT_FP")){fp++;}else{unknown++;}}
			require(p.parseInt(4)==ordered[s].size() && p.parseInt(5)==tp && p.parseInt(6)==fp && p.parseInt(8)==unknown && p.parseInt(9)==0 && p.parseInt(3)==tp+p.parseInt(7),"Native summary must agree with call rows and one-to-one truth accounting");
			counts[s][0]=tp;counts[s][1]=fp;counts[s][2]=p.parseInt(7);counts[s][3]=unknown;
		});for(boolean seen:summarySeen){require(seen,"Missing FULL_CM scope summary");}
	}
	private void write(String out){
		RrnaResourceIO.fresh(out,".mapping.tsv",".loci.tsv",".sources.tsv",".hits.tsv",".calls.tsv",".features.tsv",".precision.tsv");
		final ByteStreamWriter map=RrnaResourceIO.writer(out+".mapping.tsv"), cat=RrnaResourceIO.writer(out+".loci.tsv"), sources=RrnaResourceIO.writer(out+".sources.tsv"), hits=RrnaResourceIO.writer(out+".hits.tsv"), cr=RrnaResourceIO.writer(out+".calls.tsv"), features=RrnaResourceIO.writer(out+".features.tsv"), pr=RrnaResourceIO.writer(out+".precision.tsv");final ByteBuilder b=new ByteBuilder();
		try{
			map.println("old_review_id\tnew_review_id\tgenome\tseqid\tstrand\tstart1\tstop1");cat.println("review_id\told_review_id\tcategory\treason\tadded_scan_sources\tadded_full_GA_hits\tadded_weak_hits");
			sources.println("new_review_id\tsource_id\tmetadata_format\tmetadata\thits\tmodel\tstage\tscanned_start1\tscanned_stop1\tprojected_rows");
			hits.println("new_review_id\tsource_id\thits_file\traw_line\tmodel\tstage\tseqid\tstart1\tstop1\tstrand\tmodel_from\tmodel_to\ttrunc\tbits\tevalue\tincluded\tcall_overlap_bp\tfull_GA_overlap");
			cr.println("call_id\tgenome\tseqid\tstrand\tstart1\tstop1\temitted_family\treview_id\tcategory\toptimistic_excluded\t"+CobalaminGenomicContext.HEADER);features.println(NcrnaReviewJoin.OVERLAP_HEADER);
			for(Locus l:loci){if(l.oldId!=null){map.print(b.clear().append(l.oldId).tab().append(l.id).tab().append(l.genome).tab().append(l.seq).tab().append(l.strand).tab().append(l.start).tab().append(l.stop).nl());}
				cat.print(b.clear().append(l.id).tab().append(l.oldId==null?"NA":l.oldId).tab().append(l.category).tab().append(l.reason).tab().append(l.sources).tab().append(l.full).tab().append(l.weak).nl());
				if(l.production!=null){final Call c=l.production;cr.print(b.clear().append(c.id).tab().append(c.genome).tab().append(c.seq).tab().append(c.strand).tab().append(c.start).tab().append(c.stop).tab().append(c.family).tab().append(l.id).tab().append(l.category).tab().append(excluded(l)).tab().append(l.neighborhood).nl());for(byte[] row:l.featureRows){features.print(b.clear().append(row).nl());}}
			}
			for(byte[] row:sourceRows){sources.print(b.clear().append(row).nl());}for(byte[] row:projected){hits.print(b.clear().append(row).nl());}
			pr.println("scope\tTP\tapparent_FP\tFN\tunscored\toptimistic_excluded_union_FPs\toptimistic_remaining_FP\tpessimistic_precision\toptimistic_precision");
			for(int s=0;s<4;s++){int exclude=0;for(Call c:ordered[s]){if(c.outcome.equals("APPARENT_FP") && c.locus!=null && excluded(c.locus)){exclude++;}}
				final int tp=counts[s][0],fp=counts[s][1];require(exclude<=fp,"Optimistic exclusions cannot exceed native FPs");
				pr.print(b.clear().append(s==3?"union":MODELS[s]).tab().append(tp).tab().append(fp).tab().append(counts[s][2]).tab().append(counts[s][3]).tab().append(exclude).tab().append(fp-exclude).tab().append(tp+fp==0?"NA":String.format(Locale.ROOT,"%.8f",tp/(double)(tp+fp))).tab().append(tp+fp-exclude==0?"NA":String.format(Locale.ROOT,"%.8f",tp/(double)(tp+fp-exclude))).nl());
			}
		}finally{RrnaResourceIO.finish(map,cat,sources,hits,cr,features,pr);}
		System.err.println("COBALAMIN_PRODUCTION_PRECISION_PASS loci="+loci.size()+" productionFP="+fps+"; no added TPs, no new truth; subtype relabeling FPs retained");
	}
	private Locus locus(LineParser1 p){final Locus l=byId.get(p.parseString(0));require(l!=null,"Unknown augmented review ID");return l;}
	private static boolean sameMetadata(Locus a,Locus b){return a.key().equals(b.key()) && a.assembly.equals(b.assembly) && a.fasta.equals(b.fasta) && a.contigLength==b.contigLength && a.left==b.left && a.right==b.right && a.contextLength==b.contextLength;}
	private static boolean excluded(Locus l){return l.category.equals(REAL) || l.category.equals(GEOMETRY);}
	private static boolean category(String c){return c.equals(REAL)||c.equals(GEOMETRY)||c.equals(CONFLICT)||c.equals(PARTIAL)||c.equals(NONE);}
	private static int model(String s){for(int i=0;i<3;i++){if(MODELS[i].equals(s)){return i;}}return -1;}
	private static boolean b12(String s){
		final String value=RrnaResourceIO.optionalField(s,"Dbxref");if(value==null){return false;}
		for(int from=0;from<value.length();){int to=value.indexOf(',',from);if(to<0){to=value.length();}
			for(String model:MODELS){if(to-from==5+model.length() && value.regionMatches(true,from,"Rfam:",0,5) && value.regionMatches(from+5,model,0,model.length())){return true;}}from=to+1;
		}return false;
	}
	private static int overlap(int a,int z,int b,int y){return Math.max(0,Math.min(z,y)-Math.max(a,b)+1);}
	private static boolean bool(LineParser1 p,int col){require(p.termEquals("true",col)||p.termEquals("false",col),"Expected explicit boolean");return p.termEquals("true",col);}
	private static void read(String path,String header,java.util.function.Consumer<LineParser1> action){final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');p.set(header.getBytes(java.nio.charset.StandardCharsets.US_ASCII));final int cols=p.terms();try{RrnaResourceIO.header(in,header,path);for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==cols,"Table schema changed: "+path);action.accept(p);}}finally{require(!in.close(),"Read failed: "+path);}}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Locus {
		Locus(LineParser1 p){id=p.parseString(0);genome=p.parseString(1);assembly=p.parseString(2);fasta=p.parseString(3);seq=p.parseString(4);strand=p.parseString(6);start=p.parseInt(7);stop=p.parseInt(8);contigLength=p.parseInt(9);left=p.parseInt(10);right=p.parseInt(11);contextLength=p.parseInt(12);require(p.termEquals("cobalamin",5) && p.termEquals("PREDICTION_UNRESOLVED",13),"Expected prediction-derived cobalamin metadata");check();}
		Locus(String i,String g,String a,String s,String strand_,int start_,int stop_,int length,int left_,int right_){id=i;genome=g;assembly=a;seq=s;strand=strand_;start=start_;stop=stop_;contigLength=length;left=left_;right=right_;contextLength=length()+left+right;fasta="NA";check();}
		static Locus pilot(LineParser1 p){final String span=p.parseString(8);final int hyphen=span.indexOf('-');require(hyphen>0 && span.indexOf(',')<0 && span.indexOf(';')<0,"Historical reprojection requires one contiguous source interval");return new Locus(p.parseString(0),p.parseString(1),p.parseString(2),p.parseString(6),p.parseString(7),Integer.parseInt(span.substring(0,hyphen)),Integer.parseInt(span.substring(hyphen+1)),p.parseInt(15),p.parseInt(18),p.parseInt(19));}
		void check(){require(!id.isEmpty() && !genome.isEmpty() && (strand.equals("+")||strand.equals("-")) && start>0 && stop>=start && stop<=contigLength && left>=0 && right>=0 && contextLength==length()+left+right && contextStart()>0 && contextStop()<=contigLength,"Invalid exact locus or scanned context bounds");}
		int length(){return stop-start+1;}int contextStart(){return start-(strand.equals("+")?left:right);}int contextStop(){return stop+(strand.equals("+")?right:left);}
		String key(){return genome+'\t'+seq+'\t'+strand+'\t'+start+'\t'+stop;}
		final String id,genome,assembly,fasta,seq,strand;final int start,stop,contigLength,left,right,contextLength;String oldId,category,reason,decision,gff,gffHash,fastaHash;byte[] context,neighborhood;boolean annotatedContained;int expectedFeatures,expectedRegions,features,regions,sources,scanMask,full,fullContained,weak;Call production;
		final ObjectSet<String> featureLines=new ObjectSet<String>(String.class);final ArrayList<byte[]> featureRows=new ArrayList<byte[]>();
	}
	private static final class Call {
		Call(LineParser1 p){id=p.parseString(2);genome=p.parseString(3);seq=p.parseString(4);family=p.parseString(5);start=p.parseInt(6);stop=p.parseInt(7);strand=p.parseString(8);model=p.parseString(9);outcome=p.parseString(10);
			require(!id.isEmpty() && CobalaminProductionPrecision.model(family)>=0 && model.startsWith("cobalamin_"+family+'_') && start>0 && stop>=start && (strand.equals("+")||strand.equals("-")) && (outcome.equals("TP")||outcome.equals("APPARENT_FP")||outcome.equals("UNSCORED")),"Malformed native call");}
		String key(){return genome+'\t'+seq+'\t'+strand+'\t'+start+'\t'+stop;}boolean same(Call c){return key().equals(c.key()) && family.equals(c.family) && model.equals(c.model);}
		final String id,genome,seq,family,strand,model,outcome;final int start,stop;Locus locus;
	}
	@SuppressWarnings("unchecked")
	private CobalaminProductionPrecision(){calls=new ObjectMap[4];ordered=new ArrayList[4];for(int i=0;i<4;i++){calls[i]=new ObjectMap<String,Call>(String.class,Call.class);ordered[i]=new ArrayList<Call>();}}
	private ArrayList<Locus> loci;private int fps;
	private final ObjectMap<String,Locus> byId=new ObjectMap<String,Locus>(String.class,Locus.class),byKey=new ObjectMap<String,Locus>(String.class,Locus.class);
	private final ObjectSet<String> sourceKeys=new ObjectSet<String>(String.class);private final ArrayList<byte[]> sourceRows=new ArrayList<byte[]>(),projected=new ArrayList<byte[]>();
	private final ObjectMap<String,Call>[] calls;private final ArrayList<Call>[] ordered;private final int[][] counts=new int[4][4];private final boolean[] summarySeen=new boolean[4];
	private static final String[] MODELS={"RF00174","RF01482","RF01689"};private static final int[] CLEN={189,149,81};private static final double[] GA={31.5,50,36.5};
	private static final String REAL="REAL-B12",GEOMETRY="ANNOTATED-GEOMETRY",CONFLICT="FEATURE-CONFLICT",PARTIAL="PARTIAL/BELOW-GA",NONE="NO-SUPPORT";
	public static final String SOURCE_HEADER="new_review_id\tsource_id\tmetadata_format\tmetadata\thits\tmodel\tstage\tencoding";
	public static final String REVIEW_HEADER="review_id\tcategory\treason\tdefault_overlap_hits\tmax_overlap_hits\tfull_ga_overlap_hits\tsame_strand_B12_annotation\tdefault_weak_overlap\tdefault_weak_with_no_MAX_overlap";
	public static final String SUMMARY_HEADER="view\tscope\tgenomes\ttruth\tcalls\tTP\tFP\tFN\tunscored\texact_duplicates\tprecision\trecall\tapparent_FP_per_genome\texact5\texact3\tboth_exact";
}
