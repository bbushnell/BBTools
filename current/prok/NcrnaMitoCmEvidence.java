package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Conserves mitochondrial pilot CM evidence, physical coordinates and annotation links.
 * Canonical sha80 verification belongs to the saved producer supplying verified=.
 * Exact coordinate groups are within model AND target universe; all raw hits remain.
 * GA-support summaries require --cut_ga scans; producer verifies the exact argv,
 * and the native loader checks that flag. The reducer never filters raw rows by score.
 * Neither a GA hit nor universal preference admits truth or establishes ancestry.
 * @author Ganyu, Raiden */
public final class NcrnaMitoCmEvidence {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoCmEvidence.class,false);final Parser parser=new Parser();
		final NcrnaMitoCmEvidence x=new NcrnaMitoCmEvidence();String verified=null;
		try{
			for(String arg:pp.args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
				if(k.equals("selection")){x.selection=v;}else if(k.equals("loci")){x.locusPath=v;}else if(k.equals("models")){x.modelPath=v;}
				else if(k.equals("runsroot")){x.root=v;}else if(k.equals("members")){x.memberPath=v;}else if(k.equals("pairs")){x.pairPath=v;}
				else if(k.equals("verified")){verified=v;}else if(k.equals("outprefix")){x.prefix=v;}
				else if(k.equals("partitions")){x.partitions=Parse.parseIntKMG(v);}else if(k.equals("membership")){x.membership=v;}
				else if(k.equals("genomicroot")){x.genomicRoot=v;}else if(k.equals("genomicpartitions")){x.genomicParts=Parse.parseIntKMG(v);}
				else if(k.equals("genomicmembership")){x.genomicMembership=v;}else if(k.equals("genomiceligibility")){x.genomicEligibility=v;}
				else if(k.equals("expectedgenomicgenomes")){x.expectedGenomicGenomes=Parse.parseIntKMG(v);}else if(k.equals("expectedgenomicbases")){x.expectedGenomicBases=Parse.parseKMG(v);}
				else if(k.equals("expectedgenomes")){x.expectedGenomes=Parse.parseIntKMG(v);}else if(k.equals("expectedloci")){x.expectedLoci=Parse.parseIntKMG(v);}
				else if(k.equals("expectedmodels")){x.expectedModels=Parse.parseIntKMG(v);}else if(!parser.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}
			}
			require(verified!=null && x.selection!=null && x.locusPath!=null && x.modelPath!=null && x.root!=null && x.memberPath!=null && x.pairPath!=null && x.prefix!=null,"Require verified/selection/loci/models/runsroot/members/pairs/outprefix");
			require(x.expectedGenomes>0 && x.expectedLoci>0 && x.expectedModels>0 && !parser.append && parser.maxReads<0,"Require positive complete population counts and no append/read truncation");
				require(x.partitions>=0 && x.partitions<=128 && (x.partitions>0)==(x.membership!=null),"Partitioned scans require 1..128 partitions and the complete verified membership table; omit both for unsplit scans");
			final boolean genomicSubset=x.genomicRoot!=null;
			require(genomicSubset==(x.genomicMembership!=null) && genomicSubset==(x.genomicEligibility!=null) && genomicSubset==(x.genomicParts>0) && genomicSubset==(x.expectedGenomicGenomes>0) && genomicSubset==(x.expectedGenomicBases>0),"A restricted genomic universe requires its root, membership, eligibility, partition count and expected genome/base counts together");
			require(!genomicSubset || (x.partitions>0 && x.genomicParts<=128),"Restricted genomic evidence must retain explicitly partitioned full candidate inputs");
			ReadWrite.setZipThreads(Math.min(4,Shared.threads()));x.verified=EukPanelGffValidation.hashes(verified);
			for(String p:new String[]{x.selection,x.locusPath,x.modelPath,x.memberPath,x.pairPath}){x.bound(p);}
			x.loadSelection();x.loadLoci();x.loadMembership();x.loadModels();x.loadUniversal();x.loadRuns();x.group();x.write();
			pp.outstream.println("MITO_CM_EVIDENCE_PASS genomes="+x.genomes.size()+" loci="+x.loci.size()+" runs="+x.runs.size()+" raw_hits="+x.hits.size()+" physical_groups="+x.groups.size()+" no_truth_admission=true");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private void loadSelection(){
		read(selection,12,"accession",p->{
			final Genome g=new Genome();g.id=p.parseString(0);g.taxid=p.parseInt(1);g.organelle=p.parseString(2);g.annotation=p.parseString(3);
			g.length=p.parseInt(5);g.pad=p.parseInt(6);g.searchLength=p.parseInt(7);
			require(g.taxid>0 && g.length>0 && g.pad>=0 && g.pad<g.length && (long)g.length+g.pad==g.searchLength,"Invalid search/original coordinate map: "+g.id);
			require(g.annotation.equals("annotated") || g.annotation.equals("no_rrna"),"Unknown annotation group: "+g.id);
			require(byGenome.put(g.id,g)==null,"Repeated selected accession: "+g.id);genomes.add(g);
		});require(genomes.size()==expectedGenomes,"Selected genome count differs from expectedgenomes");
	}
	private void loadLoci(){
		locusHeader=read(locusPath,18,"locus_id",p->{
			final Locus l=new Locus();l.id=p.parseString(0);l.genome=byGenome.get(p.parseString(1));l.length=p.parseInt(16);l.raw=p.line().clone();
			require(l.genome!=null && p.parseInt(2)==l.genome.taxid && p.termEquals(l.genome.organelle,3) && p.termEquals("EXTRACTED",7) && l.length>0,"Locus must belong to exact selected source and be extracted: "+l.id);
			l.transcript=parseSegments(p.parseString(14),l.genome.length);int length=0;for(Span s:l.transcript){length+=s.hi-s.lo+1;}
			final boolean partial=NcrnaMitoPartitions.bit(p,8),pseudo=NcrnaMitoPartitions.bit(p,9);if(!partial && !pseudo){l.genome.usable++;}
			require(length==l.length,"Transcript segments must conserve extracted sequence length: "+l.id);l.footprint=normalize(l.transcript);
			require(byLocus.put(l.id,l)==null,"Duplicate annotation ID: "+l.id);loci.add(l);l.genome.loci.add(l);
		});require(loci.size()==expectedLoci,"Annotation count differs from expectedloci");
		for(Genome g:genomes){require(g.annotation.equals("no_rrna")==g.loci.isEmpty(),"Selected annotation group and retained loci disagree: "+g.id);}
	}
	private void loadModels(){
		read(modelPath,3,"accession",p->{
			final Model m=new Model();m.id=p.parseString(0);m.path=p.parseString(1);m.hash=p.parseString(2);require(bound(m.path).equals(m.hash),"Model bytes differ from verified roster");
			m.family=family(m.id);final ByteFile in=RrnaResourceIO.open(m.path);boolean acc=false;int headers=0;
			try{for(byte[] line;(line=in.nextLine())!=null;){final String s=text(line).trim();if(s.startsWith("INFERNAL1/")){headers++;}
				if(s.startsWith("ACC ")){require(s.substring(4).trim().equals(m.id),"CM accession differs from roster");acc=true;}
				else if(s.startsWith("CLEN ")){m.clen=Integer.parseInt(s.substring(5).trim());}
				else if(s.startsWith("GA ")){m.ga=Double.parseDouble(s.substring(3).trim());}
			}}finally{require(!in.close(),"CM metadata read failed");}
			require(headers==1 && acc && m.clen>0 && Double.isFinite(m.ga),"Require one complete CM and declared CLEN/GA");
			require(byModel.put(m.id,m)==null,"Duplicate model accession");models.add(m);
		});require(models.size()==expectedModels,"Model roster size mismatch");
	}
	/** Validates the partitioner's whole-accession contract before any hit is consumed. */
	private void loadMembership(){
		if(partitions==0){return;}bound(membership);
		if(genomicRoot!=null){loadGenomicEligibility();}
		final String header=read(membership,6,"universe",p->{
			final boolean genomic=p.termEquals("genomic",0);require(genomic || p.termEquals("candidate",0),"Unknown partition universe");
			final int part=p.parseInt(1);require(part>=0 && part<partitions,"Partition outside declared complete roster");
			final Genome g=byGenome.get(p.parseString(2));require(g!=null,"Partition membership names an unselected genome");
			if(genomic){
				if(genomicRoot!=null){require(p.termEquals(g.id,3) && p.parseInt(4)==g.length && p.parseInt(5)==g.searchLength && !g.originalGenomicSeen,"Original complete membership must remain unique and length-bound even when genomic searching is restricted");g.originalGenomicSeen=true;return;}
				require(p.termEquals(g.id,3) && p.parseInt(4)==g.length && p.parseInt(5)==g.searchLength && g.genomicPart<0,"Every genomic target must occur exactly once with its original and padded lengths");g.genomicPart=part;
			}else{
				final Locus l=byLocus.get(p.parseString(3));require(l!=null && l.genome==g && l.partition<0 && p.parseInt(4)==l.length && p.parseInt(5)==l.length,"Every candidate target must occur exactly once with its source accession and length");
				require(g.candidatePart<0 || g.candidatePart==part,"All annotations from one accession must stay in one candidate partition");g.candidatePart=part;l.partition=part;
			}
		});
		require(header.equals("universe\tpartition\taccession\ttarget\toriginal_length\tsearch_length"),"Require the exact native partition membership schema");
		if(genomicRoot!=null){
			bound(genomicMembership);
			final String genomicHeader=read(genomicMembership,6,"universe",p->{
				require(p.termEquals("genomic",0),"Restricted genomic membership cannot replace candidate membership");
				final Genome g=byGenome.get(p.parseString(2));final int part=p.parseInt(1);
				require(g!=null && g.genomicSelected && g.genomicPart<0 && part>=0 && part<genomicParts,"Restricted genomic target must be eligible, unique and in a declared partition");
				require(p.termEquals(g.id,3) && p.parseInt(4)==g.length && p.parseInt(5)==g.searchLength,"Restricted membership retains original whole-genome geometry");g.genomicPart=part;
			});require(genomicHeader.equals(header),"Restricted and original membership schemas must match");
		}
		for(Genome g:genomes){require((genomicRoot==null || g.originalGenomicSeen) && (g.genomicSelected?(g.genomicPart>=0):(g.genomicPart<0)),"Genomic membership must cover exactly the declared searched universe: "+g.id);}
		for(Locus l:loci){require(l.partition>=0,"Missing candidate partition membership: "+l.id);}
	}
	/** Re-derives eligibility from full source loci, so unsearched genomes cannot masquerade as negative searches. */
	private void loadGenomicEligibility(){
		assert(genomicRoot!=null && !genomes.isEmpty()):"Restricted genomic coverage is defined against the loaded full census";
		bound(genomicEligibility);
		final String header=read(genomicEligibility,7,"accession",p->{
			final Genome g=byGenome.get(p.parseString(0));require(g!=null && !g.eligibilitySeen,"Eligibility must name every source genome exactly once");g.eligibilitySeen=true;
			g.genomicSelected=NcrnaMitoPartitions.bit(p,3);
			require(p.parseInt(1)==g.loci.size() && p.parseInt(2)==g.usable && g.genomicSelected==(g.usable<2),"Genomic eligibility is fewer than two non-partial/non-pseudo source rRNA annotations, independent of product names");
			require(p.parseInt(4)==g.length && p.parseInt(5)==g.pad && p.parseInt(6)==g.searchLength,"Eligibility cannot alter the original padded genomic coordinate system");
		});require(header.equals("accession\trrna_features\tusable_rrna\tgenomic_selected\toriginal_length\tprefix_added\tsearch_length"),"Require the exact annotation-first eligibility schema");
		int count=0;long bases=0;for(Genome g:genomes){require(g.eligibilitySeen,"Missing genomic eligibility record: "+g.id);if(g.genomicSelected){count++;bases+=g.searchLength;}}
		require(count==expectedGenomicGenomes && bases==expectedGenomicBases,"Selected genomic counts and expanded bases must match the reviewed census before grading");
	}
	private void loadUniversal(){
		read(memberPath,3,"locus_id",p->{final Locus l=byLocus.get(p.parseString(0));if(l==null){return;}
			require(l.universal==null && p.parseInt(2)==l.length,"Universal member must bind once to exact annotated sequence length: "+l.id);
			final String id=p.parseString(1);Universal u=universals.get(id);if(u==null){u=new Universal(id,l.length);universals.put(id,u);}require(u.length==l.length,"Shared unique ID has inconsistent source lengths");l.universal=u;
		});for(Locus l:loci){require(l.universal!=null,"Annotation lacks a universal-comparison member: "+l.id);}
		pairHeader=read(pairPath,17,"unique_id",p->{final Universal u=universals.get(p.parseString(0));if(u==null){return;}
			final int model=p.termEquals("universal",1)?0:p.termEquals("euk18S_universal",1)?1:-1;require(model>=0 && u.rows[model]==null,"Invalid/duplicate universal pair");
			require(p.parseInt(2)==u.length,"Universal pair length differs from source annotation");u.rows[model]=p.line().clone();
			u.identities[model]=p.parseDouble(5);require(Double.isFinite(u.identities[model]) && u.identities[model]>=0 && u.identities[model]<=1,"Universal returned identity must be a finite fraction");
		});for(Locus l:loci){require(l.universal.rows[0]!=null && l.universal.rows[1]!=null,"Every annotation needs both universal comparisons");}
	}
	private void loadRuns(){
		for(Model m:models){for(String kind:new String[]{"genomic","candidate"}){for(int part=0;part<partitionCount(kind);part++){
			final Run r=new Run();r.model=m;r.kind=kind;r.partition=partitions==0?-1:part;r.id=m.id+"_"+kind+(partitions==0?"":"_p"+part);r.folder=(kind.equals("genomic") && genomicRoot!=null?genomicRoot:root)+"/"+r.id;
			r.table=r.folder+"/hits.tbl";r.hash=bound(r.table);
			for(String file:new String[]{"MANIFEST.sha80","exit_code.txt","COMPLETE.sha80","command.txt","time.txt","host.txt","INPUTS.sha80","stdout","stderr","producer.sh"}){bound(r.folder+"/"+file);}
			require(single(r.folder+"/exit_code.txt").equals("0"),"Cannot consume unfinished/failed scan: "+r.id);
			require(single(r.folder+"/COMPLETE.sha80").equals(bound(r.folder+"/MANIFEST.sha80")),"Scan seal differs from verified bytes");
			require(Arrays.asList(single(r.folder+"/command.txt").trim().split("\\s+")).contains("--cut_ga"),"GA-support summaries require a declared --cut_ga scan; exact argv is verified by the producer");
			final String jobPath=r.folder+"/job_id.txt",headPath=r.folder+"/execution.txt";
			r.execution=execution(verified.get(jobPath)==null?null:single(jobPath),verified.get(headPath)==null?null:single(headPath),single(r.folder+"/host.txt"));
			r.job=r.execution.job;
			final String[] time=single(r.folder+"/time.txt").split(" ");require(time.length==3,"Expected saved wall/user/system fields");
			r.wall=number(time[0],"wall_seconds=");r.user=number(time[1],"user_seconds=");r.system=number(time[2],"system_seconds=");
			parseTable(r.table,m.id,m.clen,m.ga,h->{h.run=r;map(h);r.hits.add(h);hits.add(h);});runs.add(r);
		}}}
		require(runs.size()==expectedModels*(partitionCount("genomic")+partitionCount("candidate")),"Every model, universe and partition must have a completed scan, including zero-hit runs");
	}
	private int partitionCount(String kind){return kind.equals("genomic") && genomicRoot!=null?genomicParts:Math.max(1,partitions);}
	/** Reuses the native 17-field parser and EukPanelCmEvidence completion/format contract. */
	public static void parseTable(String path,String model,int clen,double ga,java.util.function.Consumer<Hit> consumer){
		require(clen>0 && Double.isFinite(ga),"Pinned CM metadata required");
		final ByteFile in=RrnaResourceIO.open(path);final RrnaCmTruthBuilder.Fields f=new RrnaCmTruthBuilder.Fields();
		boolean program=false,version=false,done=false;long lineNo=0;
		try{for(byte[] line;(line=in.nextLine())!=null;){lineNo++;if(line.length==0){continue;}require(!done,"Content follows Infernal completion footer");
			if(line[0]=='#'){final String s=text(line);if(s.startsWith("# Program:")){require(s.substring(10).trim().equals("cmsearch"),"Expected cmsearch");program=true;}
				if(s.startsWith("# Version:")){require(s.substring(10).trim().startsWith("1.1.5 "),"Expected pinned Infernal1.1.5");version=true;}
				if(s.equals("# [ok]")){done=true;}continue;}
			f.set(line);require(f.eq(3,model) && f.eq(4,"cm"),"Undeclared model or HMM-only row in CM evidence");
			final Hit h=new Hit();h.line=lineNo;h.raw=line.clone();h.target=f.text(0);h.mf=f.number(5);h.mt=f.number(6);h.from=f.number(7);h.to=f.number(8);
			final String strand=f.text(9);require(strand.equals("+") || strand.equals("-"),"Invalid Infernal strand");h.strand=strand.charAt(0);
			h.trunc=f.text(10);h.inclusion=f.text(16);h.score=Parse.parseDouble(f.bytes,f.start[14],f.end[14]);h.evalue=f.text(15);
			final double evalue=Double.parseDouble(h.evalue);
			require(h.mf>=1 && h.mt>=h.mf && h.mt<=clen && h.from>0 && h.to>0,"Raw CM coordinates outside valid ranges");
			require(h.strand=='+' && h.from<=h.to || h.strand=='-' && h.from>=h.to,"Strand and target-coordinate order disagree");
			require(h.trunc.equals("no") || h.trunc.equals("5'") || h.trunc.equals("3'") || h.trunc.equals("5'&3'"),"Unknown truncation flag");
			require(h.inclusion.equals("!") || h.inclusion.equals("?"),"Unknown inclusion flag");
			require(Double.isFinite(h.score) && Double.isFinite(evalue) && evalue>=0,"Invalid numeric CM evidence");consumer.accept(h);
		}}finally{require(!in.close(),"CM table read failure: "+path);}
		require(program && version && done,"Incomplete or unrecognized Infernal table: "+path);
	}
	private void map(Hit h){
		if(h.run.kind.equals("genomic")){
			h.genome=byGenome.get(h.target);require(h.genome!=null,"Genomic hit names unselected accession: "+h.target);
			require(h.genome.genomicSelected && (partitions==0 || h.genome.genomicPart==h.run.partition),"Genomic hit belongs to an excluded genome or another partition: "+h.target);
			final NcrnaMitoPilot.MappedHit m=NcrnaMitoPilot.mapHit(h.from,h.to,h.strand,h.genome.length,h.genome.pad);
			h.appended=m.appendedDuplicate;h.crosses=m.crossesOrigin;h.spans=genomicSpans(h.from,h.to,h.strand,h.genome.length,h.genome.pad);
		}else{
			h.source=byLocus.get(h.target);require(h.source!=null,"Candidate hit names unknown annotation: "+h.target);h.genome=h.source.genome;
			require(partitions==0 || h.source.partition==h.run.partition,"Candidate hit belongs to another partition: "+h.target);
			h.spans=project(h.source.transcript,h.from,h.to,h.strand);h.crosses=h.genome.pad>0 && crosses(h.spans);
		}
		h.footprint=normalize(h.spans);h.physical=geometry(h.footprint);
		h.physicalGroupKey=physicalKey(h.run.kind,h.run.model.id,h.genome.id,h.footprint);
	}
	private void group(){
		for(Hit h:hits){
			final String key=h.physicalGroupKey;Group g=byGroup.get(key);
			if(g==null){g=new Group();g.id=groups.size();g.representative=h;groups.add(g);h.genome.physicalGroups.add(g);byGroup.put(key,g);h.run.groups++;}
			g.hits.add(h);h.group=g;
			g.representative=chooseRepresentative(g.representative,h);
		}
	}
	/** A lone appended hit is retained. Prefer an original representation only when
	 * another hit has the exact same model/universe/accession/strand/physical footprint. */
	public static Hit chooseRepresentative(Hit previous,Hit candidate){
		require(candidate!=null && candidate.physicalGroupKey!=null && Double.isFinite(candidate.score),"A retained raw hit needs a bound physical key and finite score");
		if(previous==null){return candidate;}
		require(candidate.physicalGroupKey.equals(previous.physicalGroupKey),"Cannot arbitrate different models, target universes or physical loci by score");
		if(candidate.appended!=previous.appended){return candidate.appended?previous:candidate;}
		if(candidate.score!=previous.score){return candidate.score>previous.score?candidate:previous;}
		if(candidate.mt-candidate.mf!=previous.mt-previous.mf){return candidate.mt-candidate.mf>previous.mt-previous.mf?candidate:previous;}
		return candidate.line<previous.line?candidate:previous;
	}
	private void write(){
		RrnaResourceIO.fresh(prefix,".hits.tsv.gz",".raw_rows.tsv.gz",".associations.tsv.gz",".overlaps.tsv.gz",".loci.tsv",".genomes.tsv",".runs.tsv",".execution.tsv",".universal_pairs.tsv",".cross.tsv",".cost.tsv");
		final ArrayList<ByteStreamWriter> opened=new ArrayList<ByteStreamWriter>();final ByteBuilder b=new ByteBuilder();
		try{
			final ByteStreamWriter raw=out(".raw_rows.tsv.gz",opened),hw=out(".hits.tsv.gz",opened),aw=out(".associations.tsv.gz",opened),cw=out(".overlaps.tsv.gz",opened);
			final ByteStreamWriter lw=out(".loci.tsv",opened),gw=out(".genomes.tsv",opened),rw=out(".runs.tsv",opened),uw=out(".universal_pairs.tsv",opened);
			raw.println("run_id\ttable_line\toriginal_tblout_row");
			hw.println("run_id\ttable_sha80\ttable_line\tmodel\tfamily\tuniverse\ttarget\taccession\tannotation_source\tmodel_from\tmodel_to\tclen\tmodel_coverage\traw_from\traw_to\traw_strand\ttranscript_ordered_segments\tphysical_segments\torigin_crossing\tappended_representation\texact_physical_group\trepresentative\ttrunc\tinclusion\tbits\tevalue\tga\tga_margin");
			aw.println("run_id\ttable_line\tphysical_group\trepresentative\tlocus_id\tsame_strand_overlap_bp\topposite_strand_overlap_bp\tsource_length\thit_physical_length");
			cw.println("accession\tgroup_a\tmodel_a\tuniverse_a\tgroup_b\tmodel_b\tuniverse_b\tsame_strand_overlap_bp\topposite_strand_overlap_bp\trelation");
			for(Hit h:hits){final Model m=h.run.model;final boolean representative=h.group.representative==h;
				raw.print(b.clear().append(h.run.id).tab().append(h.line).tab().append(h.raw).nl());
				hw.print(b.clear().append(h.run.id).tab().append(h.run.hash).tab().append(h.line).tab().append(m.id).tab().append(m.family).tab().append(h.run.kind)
					.tab().append(h.target).tab().append(h.genome.id).tab().append(h.source==null?"NA":h.source.id).tab().append(h.mf).tab().append(h.mt).tab().append(m.clen)
					.tab().append((double)(h.mt-h.mf+1)/m.clen,9).tab().append(h.from).tab().append(h.to).tab().append(h.strand).tab().append(geometry(h.spans)).tab().append(h.physical)
					.tab().append(h.crosses).tab().append(h.appended).tab().append(h.group.id).tab().append(representative).tab().append(h.trunc).tab().append(h.inclusion)
					.tab().append(h.score,9).tab().append(h.evalue).tab().append(m.ga,9).tab().append(h.score-m.ga,9).nl());
				if(representative){h.genome.groups++;}
				for(Locus l:h.genome.loci){final long same=overlap(h.footprint,l.footprint,true),opposite=overlap(h.footprint,l.footprint,false);
					if(same+opposite==0){continue;}
					aw.print(b.clear().append(h.run.id).tab().append(h.line).tab().append(h.group.id).tab().append(representative).tab().append(l.id).tab().append(same).tab().append(opposite).tab().append(l.length).tab().append(spanLength(h.footprint)).nl());
					if(representative){if(same>0){l.support[h.run.kind.equals("genomic")?0:1][m.family.equals("SSU")?0:1]++;}if(opposite>0){l.opposite++;}}
				}
			}
			// Per-genome lists preserve global insertion order, so this skips unrelated
			// accessions without changing the old global (i,j) output order.
			for(Group left:groups){final Hit a=left.representative;for(Group right:a.genome.physicalGroups){
				if(right.id<=left.id){continue;}final Hit z=right.representative;
				assert(a.genome==z.genome):"Physical overlap candidates must share a genome; cross-accession coordinates are unrelated";
				final long same=overlap(a.footprint,z.footprint,true),opposite=overlap(a.footprint,z.footprint,false);if(same+opposite==0){continue;}
				final String relation=a.run.model.family.equals(z.run.model.family)?"SAME_FAMILY_EVIDENCE":"SSU_LSU_OVERLAP_REVIEW";
				cw.print(b.clear().append(a.genome.id).tab().append(a.group.id).tab().append(a.run.model.id).tab().append(a.run.kind).tab().append(z.group.id).tab().append(z.run.model.id).tab().append(z.run.kind).tab().append(same).tab().append(opposite).tab().append(relation).nl());
			}}
			lw.println("locus_id\taccession\tunique_id\tgenomic_ssu_groups\tgenomic_lsu_groups\tcandidate_ssu_groups\tcandidate_lsu_groups\topposite_strand_groups\tstatus\tsource_"+locusHeader.replace("\t","\tsource_")+(genomicRoot==null?"":"\tgenomic_search"));
			uw.println("locus_id\t"+pairHeader);
			for(Locus l:loci){final int ssu=l.support[0][0]+l.support[1][0],lsu=l.support[0][1]+l.support[1][1];
				final String status=ssu==0 && lsu==0?"NO_SAME_STRAND_GA_HIT_UNRESOLVED":ssu>0 && lsu>0?"SSU_LSU_CONFLICT_REVIEW":ssu>0?"SSU_SUPPORT_REVIEW":"LSU_SUPPORT_REVIEW";
				b.clear().append(l.id).tab().append(l.genome.id).tab().append(l.universal.id).tab();
				if(l.genome.genomicSelected){b.append(l.support[0][0]).tab().append(l.support[0][1]);}else{b.append("NA\tNA");}
				b.tab().append(l.support[1][0]).tab().append(l.support[1][1]).tab().append(l.opposite).tab().append(status).tab().append(l.raw);
				if(genomicRoot!=null){b.tab().append(l.genome.genomicSelected?"SEARCHED":"NOT_SEARCHED");}lw.print(b.nl());
				for(byte[] row:l.universal.rows){uw.print(b.clear().append(l.id).tab().append(row).nl());}
			}
			gw.println("accession\ttaxid\torganelle_class\tannotation_group\tannotation_loci\tphysical_groups_across_separate_universes\tstatus"+(genomicRoot==null?"":"\tgenomic_search\tcandidate_search"));
			for(Genome g:genomes){b.clear().append(g.id).tab().append(g.taxid).tab().append(g.organelle).tab().append(g.annotation).tab().append(g.loci.size()).tab().append(g.groups).tab().append(g.groups==0?"NO_GA_HIT_UNRESOLVED":g.loci.isEmpty()?"CM_ONLY_REVIEW":"ANNOTATED_REVIEW");
				if(genomicRoot!=null){b.tab().append(g.genomicSelected?"SEARCHED":"NOT_SEARCHED").tab().append(g.loci.isEmpty()?"NO_ANNOTATED_TARGET":"SEARCHED");}gw.print(b.nl());}
			rw.println("run_id\tjob_id\tmodel\tfamily\tuniverse\ttable\ttable_sha80\tclen\tga\traw_hits\texact_physical_groups\twall_seconds\tuser_seconds\tsystem_seconds");
			for(Run r:runs){rw.print(b.clear().append(r.id).tab().append(r.job).tab().append(r.model.id).tab().append(r.model.family).tab().append(r.kind).tab().append(r.table).tab().append(r.hash).tab().append(r.model.clen).tab().append(r.model.ga,3).tab().append(r.hits.size()).tab().append(r.groups).tab().append(r.wall,3).tab().append(r.user,3).tab().append(r.system,3).nl());}
			final ByteStreamWriter ew=out(".execution.tsv",opened);ew.println("run_id\texecutor\tscheduler_job_id\thost\tprocess_id");
			for(Run r:runs){final Execution e=r.execution;ew.print(b.clear().append(r.id).tab().append(e.executor).tab().append(e.job).tab().append(e.host).tab().append(e.pid).nl());}
		}finally{RrnaResourceIO.finish(opened.toArray(new ByteStreamWriter[opened.size()]));}
		writeSummary();
	}
	/** Annotation-weighted evidence/preference cross-table, not biological subtype truth.
	 * Sum of job wall times is explicitly distinct from parallel elapsed wall time. */
	private void writeSummary(){
		final int scopes=genomicRoot==null?1:2;final int[][][] counts=new int[scopes][4][3];final double[][][][] identities=new double[scopes][4][3][2];
		final int[] denominators=new int[scopes];int total=0;
		for(Locus l:loci){final int ssu=l.support[0][0]+l.support[1][0],lsu=l.support[0][1]+l.support[1][1];
			final int status=ssu==0 && lsu==0?0:ssu>0 && lsu>0?3:ssu>0?1:2;
			final double a=l.universal.identities[0],b=l.universal.identities[1];final int preference=a>b?0:a==b?1:2;
			final int scope=l.genome.genomicSelected?0:1;
			counts[scope][status][preference]++;identities[scope][status][preference][0]+=a;identities[scope][status][preference][1]+=b;denominators[scope]++;total++;
		}
		require(total==expectedLoci,"Cross-table must conserve every original annotation, including unresolved cases");
		final ByteStreamWriter cross=RrnaResourceIO.writer(prefix+".cross.tsv"),cost=RrnaResourceIO.writer(prefix+".cost.tsv");final ByteBuilder b=new ByteBuilder();
		try{
			cross.println((genomicRoot==null?"":"genomic_search\t")+"cm_annotation_status\tuniversal_preference\tannotations\tmean_quantum_16S\tmean_quantum_18S\tmean_18S_minus_16S\tannotation_denominator");
			for(int scope=0;scope<scopes;scope++){for(int s=0;s<4;s++){for(int p=0;p<3;p++){final int n=counts[scope][s][p];b.clear();
				if(genomicRoot!=null){b.append(scope==0?"SEARCHED":"NOT_SEARCHED").tab();}b.append(CM_STATUSES[s]).tab().append(PREFERENCES[p]).tab().append(n).tab();
				if(n==0){b.append("NA\tNA\tNA");}else{final double a=identities[scope][s][p][0]/n,c=identities[scope][s][p][1]/n;b.append(a,9).tab().append(c,9).tab().append(c-a,9);}
				cross.print(b.tab().append(denominators[scope]).nl());
			}}}
			cost.println(partitions==0?"universe\ttarget_sequences\ttarget_bases\truns\tsum_job_wall_seconds\tmax_job_wall_seconds\tsum_user_seconds\tsum_system_seconds\tsum_cpu_seconds\tmin_job_wall_seconds":"universe\ttarget_sequences\ttarget_bases\truns\tsum_task_wall_seconds\tmax_task_wall_seconds\tsum_user_seconds\tsum_system_seconds\tsum_cpu_seconds\tmin_task_wall_seconds");
			for(String kind:new String[]{"genomic","candidate"}){int n=0;double wall=0,max=0,user=0,system=0,min=Double.MAX_VALUE;long bases=0;
				for(Run r:runs){if(r.kind.equals(kind)){n++;wall+=r.wall;max=Math.max(max,r.wall);min=Math.min(min,r.wall);user+=r.user;system+=r.system;}}
				int targets=0;
				if(kind.equals("genomic")){for(Genome g:genomes){if(g.genomicSelected){bases+=g.searchLength;targets++;}}}else{for(Locus l:loci){bases+=l.length;targets++;}}
				require(n==expectedModels*partitionCount(kind) && bases>0,"Every cost row needs the full model/partition roster and exact searched target universe counted once");
				cost.print(b.clear().append(kind).tab().append(targets).tab().append(bases).tab().append(n).tab().append(wall,3).tab().append(max,3).tab().append(user,3).tab().append(system,3).tab().append(user+system,3).tab().append(min,3).nl());
			}
		}finally{RrnaResourceIO.finish(cross,cost);}
	}
	public static ArrayList<Span> parseSegments(String value,int genomeLength){
		require(value!=null && !value.isEmpty() && genomeLength>0,"Transcript geometry and source length required");final ArrayList<Span> result=new ArrayList<Span>();
		for(String token:value.split(",")){final int dots=token.indexOf(".."),colon=token.indexOf(':');require(dots>0 && colon>dots+2 && colon==token.length()-2,"Malformed retained transcript segment: "+token);
			final int lo=Integer.parseInt(token.substring(0,dots)),hi=Integer.parseInt(token.substring(dots+2,colon));final char strand=token.charAt(colon+1);
			require(lo>=1 && hi>=lo && hi<=genomeLength && (strand=='+' || strand=='-'),"Transcript segment out of range");result.add(new Span(lo,hi,strand));
		}return result;
	}
	/** Projects a hit on a concatenated, sense-oriented candidate through every original segment. */
	public static ArrayList<Span> project(ArrayList<Span> transcript,int from,int to,char strand){
		require(transcript!=null && !transcript.isEmpty() && from>0 && to>0 && (strand=='+' && from<=to || strand=='-' && from>=to),"Invalid candidate projection input");
		final int lo=Math.min(from,to),hi=Math.max(from,to);int offset=1;final ArrayList<Span> result=new ArrayList<Span>();
		for(Span s:transcript){final int end=offset+s.hi-s.lo;final int a=Math.max(lo,offset),b=Math.min(hi,end);
			if(a<=b){final int x=s.strand=='+'?s.lo+a-offset:s.hi-(a-offset),y=s.strand=='+'?s.lo+b-offset:s.hi-(b-offset);
				result.add(new Span(Math.min(x,y),Math.max(x,y),strand=='+'?s.strand:s.strand=='+'?'-':'+'));}
			offset=end+1;
		}require(hi<offset && spanLength(result)==(long)hi-lo+1,"Projected segments must conserve the complete CM target interval");
		if(strand=='-'){Collections.reverse(result);}return result;
	}
	public static ArrayList<Span> genomicSpans(int from,int to,char strand,int length,int pad){
		final NcrnaMitoPilot.MappedHit m=NcrnaMitoPilot.mapHit(from,to,strand,length,pad);final ArrayList<Span> spans=new ArrayList<Span>();
		if(!m.crossesOrigin){spans.add(new Span(Math.min(m.from,m.to),Math.max(m.from,m.to),strand));}
		else{final int lo=Math.min(from,to),hi=Math.max(from,to);spans.add(new Span(lo,length,strand));spans.add(new Span(1,hi-length,strand));if(strand=='-'){Collections.reverse(spans);}}
		require(spanLength(spans)==(long)Math.abs(from-to)+1,"Circular projection must conserve hit span");return spans;
	}
	public static ArrayList<Span> normalize(ArrayList<Span> source){
		require(source!=null && !source.isEmpty(),"Physical footprint cannot be empty");final ArrayList<Span> sorted=new ArrayList<Span>(source),result=new ArrayList<Span>();
		sorted.sort(Comparator.comparingInt((Span s)->s.strand).thenComparingInt(s->s.lo).thenComparingInt(s->s.hi));
		for(Span s:sorted){if(result.isEmpty()){result.add(s);continue;}final Span old=result.get(result.size()-1);
			if(old.strand==s.strand && (long)s.lo<=old.hi+1L){result.set(result.size()-1,new Span(old.lo,Math.max(old.hi,s.hi),s.strand));}else{result.add(s);}}
		return result;
	}
	public static long overlap(ArrayList<Span> a,ArrayList<Span> b,boolean sameStrand){
		require(a!=null && b!=null,"Both normalized footprints are required");long count=0;
		for(Span x:a){for(Span y:b){if((x.strand==y.strand)==sameStrand){count+=Math.max(0,Math.min(x.hi,y.hi)-Math.max(x.lo,y.lo)+1);}}}return count;
	}
	public static String geometry(ArrayList<Span> spans){final ByteBuilder b=new ByteBuilder();for(Span s:spans){if(b.length()>0){b.append(',');}b.append(s.lo).append("..").append(s.hi).append(':').append(s.strand);}return b.toString();}
	public static String physicalKey(String universe,String model,String accession,ArrayList<Span> spans){require(universe!=null && model!=null && accession!=null,"Physical groups require explicit target universe, model and accession");return universe+"\t"+model+"\t"+accession+"\t"+geometry(normalize(spans));}
	private static boolean crosses(ArrayList<Span> spans){for(int i=1;i<spans.size();i++){final Span a=spans.get(i-1),b=spans.get(i);if(a.strand==b.strand && (a.strand=='+' && b.lo<a.lo || a.strand=='-' && b.hi>a.hi)){return true;}}return false;}
	/** Both inputs are optional individually, but exactly one verified execution source is required.
	 * A head PID is never put in the scheduler-job column; original Slurm rows remain unchanged. */
	public static Execution execution(String job,String head,String host){
		require((job==null)!=(head==null),"Exactly one bound scheduler job or head execution record is required");
		require(host!=null && host.matches("[A-Za-z0-9][A-Za-z0-9.-]*"),"Recorded execution hostname required");
		if(job!=null){require(job.matches("[1-9][0-9]*"),"Invalid scheduler job ID");return new Execution("slurm",job,host,"NA");}
		require(head.matches("executor=head pid=[1-9][0-9]*"),"Head execution must explicitly identify a positive process ID");
		return new Execution("head","NA",host,head.substring("executor=head pid=".length()));
	}
	private static long spanLength(ArrayList<Span> spans){long n=0;for(Span s:spans){n+=(long)s.hi-s.lo+1;}return n;}
	private String bound(String path){final String hash=verified.get(path);require(hash!=null,"Input missing from canonical-sha80-verified manifest: "+path);return hash;}
	private static String single(String path){final ByteFile in=RrnaResourceIO.open(path);try{final byte[] line=in.nextLine();require(line!=null && in.nextLine()==null,"Expected one line: "+path);return text(line);}finally{require(!in.close(),"Read failed: "+path);}}
	private static double number(String text,String prefix){require(text.startsWith(prefix),"Unexpected timing field");final double n=Double.parseDouble(text.substring(prefix.length()));require(Double.isFinite(n) && n>=0,"Invalid observed time");return n;}
	private static String family(String id){if(id.equals("RF00177") || id.equals("RF01959") || id.equals("RF01960") || id.equals("RF02542") || id.equals("RF02545")){return "SSU";}if(id.equals("RF02540") || id.equals("RF02541") || id.equals("RF02543") || id.equals("RF02546")){return "LSU";}throw new IllegalArgumentException("Undeclared SSU/LSU model: "+id);}
	private static String text(byte[] line){return new String(line,StandardCharsets.UTF_8);}
	private static String read(String path,int columns,String first,java.util.function.Consumer<LineParser1> row){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final String observed;
		try{final byte[] header=in.nextLine();require(header!=null,"Missing table header: "+path);observed=text(header);p.set(header);require(p.terms()==columns && p.termEquals(first,0),"Unexpected schema: "+path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==columns,"Unexpected row width: "+path);row.accept(p);}}
		finally{require(!in.close(),"Table read failed: "+path);}
		return observed;
	}
	private ByteStreamWriter out(String suffix,ArrayList<ByteStreamWriter> opened){final ByteStreamWriter w=RrnaResourceIO.writer(prefix+suffix);opened.add(w);return w;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final class Span {public Span(int l,int h,char s){require(l>0 && h>=l && (s=='+' || s=='-'),"Invalid physical segment");lo=l;hi=h;strand=s;}public final int lo,hi;public final char strand;}
	public static final class Hit {public int mf,mt,from,to;public char strand;public long line;public double score;public String target,trunc,inclusion,evalue,physicalGroupKey;public boolean appended;byte[] raw;Run run;Genome genome;Locus source;Group group;ArrayList<Span> spans,footprint;String physical;boolean crosses;}
	private static final class Genome {String id,organelle,annotation;int taxid,length,pad,searchLength,groups,usable;int genomicPart=-1,candidatePart=-1;boolean genomicSelected=true,eligibilitySeen,originalGenomicSeen;final ArrayList<Locus> loci=new ArrayList<Locus>();final ArrayList<Group> physicalGroups=new ArrayList<Group>();}
	private static final class Locus {String id;Genome genome;int length,opposite;int partition=-1;byte[] raw;ArrayList<Span> transcript,footprint;Universal universal;final int[][] support=new int[2][2];}
	private static final class Universal {Universal(String s,int n){id=s;length=n;}final String id;final int length;final byte[][] rows=new byte[2][];final double[] identities=new double[2];}
	private static final class Model {String id,path,hash,family;int clen;double ga=Double.NaN;}
	public static final class Execution {Execution(String e,String j,String h,String p){executor=e;job=j;host=h;pid=p;}public final String executor,job,host,pid;}
	private static final class Run {String id,kind,folder,table,hash,job;Model model;Execution execution;double wall,user,system;int groups,partition;final ArrayList<Hit> hits=new ArrayList<Hit>();}
	private static final class Group {int id;Hit representative;final ArrayList<Hit> hits=new ArrayList<Hit>();}
	private String selection,locusPath,modelPath,root,memberPath,pairPath,prefix,locusHeader,pairHeader,membership;private int expectedGenomes=16,expectedLoci=22,expectedModels=9,partitions=0;
	private String genomicRoot,genomicMembership,genomicEligibility;private int genomicParts,expectedGenomicGenomes;private long expectedGenomicBases;
	private ObjectMap<String,String> verified;private final ObjectMap<String,Genome> byGenome=new ObjectMap<String,Genome>(String.class,Genome.class);
	private final ObjectMap<String,Locus> byLocus=new ObjectMap<String,Locus>(String.class,Locus.class);private final ObjectMap<String,Model> byModel=new ObjectMap<String,Model>(String.class,Model.class);
	private final ObjectMap<String,Universal> universals=new ObjectMap<String,Universal>(String.class,Universal.class);private final ObjectMap<String,Group> byGroup=new ObjectMap<String,Group>(String.class,Group.class);
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();private final ArrayList<Locus> loci=new ArrayList<Locus>();private final ArrayList<Model> models=new ArrayList<Model>();
	private final ArrayList<Run> runs=new ArrayList<Run>();private final ArrayList<Hit> hits=new ArrayList<Hit>();private final ArrayList<Group> groups=new ArrayList<Group>();
	private static final String[] CM_STATUSES={"NO_SAME_STRAND_GA_HIT_UNRESOLVED","SSU_SUPPORT_REVIEW","LSU_SUPPORT_REVIEW","SSU_LSU_CONFLICT_REVIEW"};
	private static final String[] PREFERENCES={"16S","TIE","18S"};
}
