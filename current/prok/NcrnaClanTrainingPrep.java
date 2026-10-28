package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Routes every admitted core to a training group without changing biological labels.
 * Keeps the accepted RF00174 product intact and attributes every new admission.
 * @author Brian, Raiden
 */
public final class NcrnaClanTrainingPrep {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaClanTrainingPrep.class, false);
		String manifest=null, out=null;int expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("inputs")){manifest=v;}else if(k.equals("outprefix")){out=v;}else if(k.equals("expected")){expected=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(manifest!=null && out!=null && expected>0, "Bind all three complete admitted cohorts and the combined denominator");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{final NcrnaClanTrainingPrep x=new NcrnaClanTrainingPrep();x.inputs(manifest, expected);x.write(out, expected);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void inputs(String path, int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int total=0;
		try{
			RrnaResourceIO.header(in, "family\tfasta\toverlay\tprevious_overlay\tjoined\tselection\trawmodels\texpected", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==8 && cohorts.size()<3 && p.termEquals(IDS[cohorts.size()], 0), "Cohorts must be RF00174, RF01482 and RF01689 in canonical order");
				final Cohort c=new Cohort(p);require(c.fasta.startsWith("/") && c.expected>0, "Each complete cohort needs an absolute FASTA and positive denominator");
				cohorts.add(c);total+=c.expected;
			}
		}finally{require(!in.close(), "Cohort manifest read failed");}
		require(cohorts.size()==3 && total==expected, "The combined training population must equal all three declared cohorts");
		for(Cohort c:cohorts){
			if(c.family.equals(IDS[0])){require(c.overlay.equals("NA") && c.previous.equals("NA") && c.joined.equals("NA") && c.selection.equals("NA") && c.rawmodels.equals("NA"), "RF00174 reuses its explicitly accepted product, not a fabricated per-locus CM overlay");}
			else{overlay(c, false);overlay(c, true);joined(c);selection(c);}
		}
	}
	private void overlay(Cohort c, boolean previous){
		final ByteFile in=RrnaResourceIO.open(previous?c.previous:c.overlay);final LineParser1 p=new LineParser1('\t');int count=0, admitted=0;
		try{
			final byte[] header=in.nextLine();require(header!=null, "Overlay schema is required");
			final String text=new String(header, java.nio.charset.StandardCharsets.UTF_8), current=NcrnaClanOverlay.HEADER;
			final boolean legacy=text.equals(current.substring(0, current.lastIndexOf('\t')));
			require(text.equals(current) || (previous && legacy), "Only the exact current or explicitly versioned first-pass overlay schema is supported");
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==(legacy?32:33) && p.termEquals(c.family, 21), "Overlay family or field count changed");
				final String id=p.parseString(0);require(id.startsWith(c.family+"_"), "Foreign source ID: "+id);
				final boolean accepted=flag(p, 23);Member m=c.members.get(id);
				if(previous){
					require(m!=null && !m.previousSeen && Arrays.equals(m.source, sourcePrefix(row)), "Previous/current overlays must preserve every original source row exactly once: "+id);
					m.previousSeen=true;m.previousAccepted=accepted;
					require(!accepted || m.accepted, "The revised support rule must not silently remove an accepted core: "+id);
				}else{require(m==null, "Repeated overlay source: "+id);m=new Member(p, row);c.members.put(id, m);c.order.add(m);}
				count++;if(accepted){admitted++;}
			}
		}finally{require(!in.close(), "Overlay read failed");}
		require(previous?count==c.order.size():admitted==c.expected, "Overlay source/admitted denominator changed");
	}
	private void joined(Cohort c){
		final ByteFile in=RrnaResourceIO.open(c.joined);final LineParser1 p=new LineParser1('\t');int count=0;
		try{
			final byte[] header=in.nextLine();require(header!=null, "Joined evidence schema required");p.set(header);
			require(p.terms()==72 && p.termEquals("source_family", 21), "Expected the complete three-model evidence schema");
			for(int m=0; m<3; m++){require(p.termEquals(IDS[m]+"_cm_state", 26+15*m) && p.termEquals(IDS[m]+"_cm_bits", 34+15*m), "Joined model order changed");}
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==72, "Joined model row is incomplete");final Member m=c.members.get(p.parseString(0));
				require(m!=null && !m.joinedSeen && Arrays.equals(m.source, sourcePrefix(row)), "Joined evidence changed or repeated source provenance");m.joinedSeen=true;count++;
				if(!m.accepted){continue;}
				int supports=0, best=-1;double score=Double.NEGATIVE_INFINITY;boolean tied=false;
				for(int i=0; i<3; i++){
					final int col=26+15*i;if(!p.termEquals("CM_FULL_SPAN_ANY_HIT_SUPPORTED", col)){continue;}supports++;
					final double bits=p.parseDouble(col+8);require(Double.isFinite(bits), "Operational routing requires finite source CM scores");
					if(bits>score){score=bits;best=i;tied=false;}else if(bits==score){tied=true;}
				}
				require(supports>0, "An admitted core must retain full raw CM support");
				if(m.memberFamily.equals("B12clade")){
					require(!tied, "Equal-score cross-labeled routing requires a recorded operational tie-break: "+m.id);
					m.group=best;m.routing=supports==1?"SOLE_SUPPORTING_CM":"HIGHEST_BITS_OPERATIONAL_ONLY";m.routingBits=score;
				}else{
					m.group=index(c.family);require(p.termEquals("CM_FULL_SPAN_ANY_HIT_SUPPORTED", 26+15*m.group), "Own-subtype members need own full CM support");
					m.routing="OWN_SUBTYPE";m.routingBits=p.parseDouble(34+15*m.group);
				}
			}
		}finally{require(!in.close(), "Joined evidence read failed");}
		require(count==c.order.size(), "Joined evidence must conserve the entire original population");
	}
	private void selection(Cohort c){
		final ObjectMap<String,String> tables=new ObjectMap<String,String>(String.class, String.class);final LineParser1 p=new LineParser1('\t');
		ByteFile in=RrnaResourceIO.open(c.rawmodels);
		try{
			RrnaResourceIO.header(in, "model_id\tclen\tga_bits\tstage\thits", c.rawmodels);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==5 && index(p.parseString(0))>=0 && (p.termEquals("DEFAULT", 3) || p.termEquals("SENSITIVE", 3)) && p.termStartsWith("/", 4), "Raw table model, stage and absolute path must remain explicit");
				final String key=p.parseString(0)+'/'+p.parseString(3);require(!tables.contains(key), "Repeated raw table binding");tables.put(key, p.parseString(4));
			}
		}finally{require(!in.close(), "Raw table manifest read failed");}
		require(tables.size()==6, "Both scan stages must be retained for all three models");
		in=RrnaResourceIO.open(c.selection);int count=0;final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{
			RrnaResourceIO.header(in, "locus_id\tmodel_id\tdefault_raw\tsensitive_raw\tdefault_overlapping\tsensitive_overlapping\tdefault_qualified\tsensitive_qualified\tselected_stage\tselected_table_line\tseq_from\tseq_to\tgeometry\tstatus", c.selection);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==14, "Raw selection row changed");final String id=p.parseString(0), model=p.parseString(1);
				final Member m=c.members.get(id);require(m!=null && index(model)>=0 && seen.add(id+'/'+model), "Raw selection must retain each source/model pair exactly once");count++;
				if(!m.accepted || !m.coordinateModel.equals(model)){continue;}
				require(p.termEquals("FULL_SUPPORT", 13), "Coordinate model must point to an actual selected full hit");
				m.stage=p.parseString(8);m.rawTable=tables.get(model+'/'+m.stage);m.rawLine=p.parseInt(9);m.geometry=p.parseString(12);
				m.defaultOverlap=p.parseInt(4);m.sensitiveOverlap=p.parseInt(5);require(m.rawTable!=null && m.rawLine>0, "A gain must resolve to its actual raw table and line");
			}
		}finally{require(!in.close(), "Raw selection read failed");}
		require(count==3*c.order.size(), "Raw selection must preserve all three models for every source");
	}
	private void write(String out, int expected){
		RrnaResourceIO.fresh(out, ".clan.fa", ".members.tsv", ".gains.tsv", ".summary.tsv");for(String id:IDS){RrnaResourceIO.fresh(out, "."+id+".fa");}
		final ByteStreamWriter[] groups=new ByteStreamWriter[3];for(int i=0; i<3; i++){groups[i]=RrnaResourceIO.writer(out+'.'+IDS[i]+".fa");}
		final ByteStreamWriter clan=RrnaResourceIO.writer(out+".clan.fa"), members=RrnaResourceIO.writer(out+".members.tsv"), gains=RrnaResourceIO.writer(out+".gains.tsv");
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);final ByteBuilder b=new ByteBuilder();final int[] counts=new int[3], cross=new int[3], gainCount=new int[3];
		try{
			members.println("locus_id\tsource_family\tmember_family\tbuild_group\trouting_rule\trouting_bits\tlength\tinput_header");
			gains.println("locus_id\tsource_family\tcoordinate_model\tselected_stage\traw_table\traw_table_line\tgeometry\tdefault_overlap_events\tsensitive_overlap_events\toriginal_interval1\tcurrent_start1\tcurrent_stop1\ttruth_subtype");
			for(Cohort c:cohorts){
				final int[] observed={0};
				RrnaResourceIO.stream(c.fasta, r->{
					final String id=RrnaResourceIO.first(r.id);require(id.startsWith(c.family+"_") && seen.add(id) && r.bases!=null && r.length()>=17 && r.id.indexOf('\t')<0,
						"Training inputs must retain unique source identities, nonempty seed-length sequences and TSV-safe headers");
					for(byte base:r.bases){require(AminoAcid.baseToNumber[base&255]>=0, "Ambiguous bases cannot enter a verified core training set: "+id);}
					final Member m=c.members.get(id);final boolean old=c.family.equals(IDS[0]);
					if(old){require(RrnaResourceIO.optionalField(r.id, "member_family")==null, "Accepted RF00174 base must not already carry a rewritten biological member label");}
					require(old || m!=null && m.accepted && !m.fastaSeen && m.length==r.length() && m.rawTable!=null, "Admitted FASTA must match exactly one fully traced overlay core: "+id);
					if(!old){
						require(RrnaResourceIO.field(r.id, "member_family").equals(m.memberFamily) && RrnaResourceIO.field(r.id, "source_family").equals(c.family)
							&& RrnaResourceIO.field(r.id, "coordinate_model").equals(m.coordinateModel)
							&& Integer.parseInt(RrnaResourceIO.field(r.id, "start"))==m.start && Integer.parseInt(RrnaResourceIO.field(r.id, "stop"))==m.stop,
							"Training core header changed its reviewed biological/coordinate provenance: "+id);m.fastaSeen=true;
					}
					final int group=old?0:m.group;final String biological=old?IDS[0]:m.memberFamily, rule=old?"ACCEPTED_RF00174_PRODUCT":m.routing;
					require(group>=0 && group<3 && RrnaResourceIO.optionalField(r.id, "build_group")==null, "Operational group must be explicit and not applied twice");
					observed[0]++;counts[group]++;if(biological.equals("B12clade")){cross[group]++;}
					b.clear().append('>').append(r.id);if(old){b.append(" member_family=").append(biological);}b.append(" build_group=").append(IDS[group]).nl().append(r.bases).nl();
					groups[group].print(b);clan.print(b);
					members.print(b.clear().append(id).tab().append(c.family).tab().append(biological).tab().append(IDS[group]).tab().append(rule).tab()
						.append(old?"NA":Double.toString(m.routingBits)).tab().append(r.length()).tab().append(r.id).nl());
					if(!old && !m.previousAccepted){
						gainCount[index(c.family)]++;gains.print(b.clear().append(id).tab().append(c.family).tab().append(m.coordinateModel).tab().append(m.stage).tab().append(m.rawTable)
							.tab().append(m.rawLine).tab().append(m.geometry).tab().append(m.defaultOverlap).tab().append(m.sensitiveOverlap).tab().append(m.originalInterval)
							.tab().append(m.start).tab().append(m.stop).tab().append(m.truth).nl());
					}
				});
				require(observed[0]==c.expected, "Every admitted member must be routed: "+c.family+" observed="+observed[0]+" expected="+c.expected);
				for(Member m:c.order){require(!m.accepted || m.fastaSeen, "An admitted overlay core is missing from the input FASTA: "+m.id);}
			}
			require(seen.size()==expected && counts[0]+counts[1]+counts[2]==expected, "Three training groups must conserve the complete clan population");
		}finally{RrnaResourceIO.finish(groups[0], groups[1], groups[2], clan, members, gains);}
		final ByteStreamWriter summary=RrnaResourceIO.writer(out+".summary.tsv");
		try{summary.println("family\ttraining_members\tclan_only_members\tnew_admissions_from_source_family");for(int i=0; i<3; i++){
			summary.print(b.clear().append(IDS[i]).tab().append(counts[i]).tab().append(cross[i]).tab().append(gainCount[i]).nl());}}
		finally{RrnaResourceIO.finish(summary);}
		System.err.println("NCRNA_CLAN_TRAINING_PREP_PASS members="+expected+"; operational routing preserves biological member labels");
	}
	private static byte[] sourcePrefix(byte[] row){
		int n=0;for(int i=0; i<row.length; i++){if(row[i]=='\t' && ++n==21){return Arrays.copyOf(row, i);}}throw new IllegalArgumentException("Missing original source prefix");
	}
	private static boolean flag(LineParser1 p, int column){require(p.termEquals("true", column) || p.termEquals("false", column), "Admission flags must be explicit");return p.termEquals("true", column);}
	private static int index(String id){for(int i=0; i<3; i++){if(id.equals(IDS[i])){return i;}}return -1;}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Cohort {
		Cohort(LineParser1 p){family=p.parseString(0);fasta=p.parseString(1);overlay=p.parseString(2);previous=p.parseString(3);joined=p.parseString(4);selection=p.parseString(5);rawmodels=p.parseString(6);expected=p.parseInt(7);}
		final String family, fasta, overlay, previous, joined, selection, rawmodels;final int expected;
		final ObjectMap<String,Member> members=new ObjectMap<String,Member>(String.class, Member.class);final ArrayList<Member> order=new ArrayList<Member>();
	}
	private static final class Member {
		Member(LineParser1 p, byte[] row){
			id=p.parseString(0);source=sourcePrefix(row);accepted=flag(p, 23);truth=p.parseString(22);memberFamily=p.parseString(31);originalInterval=p.parseString(8);coordinateModel=p.parseString(25);
			if(accepted){start=p.parseInt(27);stop=p.parseInt(28);length=p.parseInt(29);require(start>0 && stop-start+1==length && (truth.equals(p.parseString(21)) && memberFamily.equals(truth)
				|| truth.equals("CROSS_LABELED") && memberFamily.equals("B12clade")), "Accepted biological labels and corrected geometry must agree");}
			else{start=stop=length=0;require(truth.equals("UNRESOLVED") && memberFamily.equals("UNRESOLVED"), "Unadmitted sources must remain explicitly unresolved");}
		}
		final String id, truth, memberFamily, originalInterval, coordinateModel;final byte[] source;final boolean accepted;final int start, stop, length;
		boolean previousSeen, previousAccepted, joinedSeen, fastaSeen;int group=-1, rawLine, defaultOverlap, sensitiveOverlap;double routingBits;
		String routing, stage, rawTable, geometry;
	}
	private final ArrayList<Cohort> cohorts=new ArrayList<Cohort>();
	private static final String[] IDS={"RF00174", "RF01482", "RF01689"};
}
