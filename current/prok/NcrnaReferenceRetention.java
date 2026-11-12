package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Comparator;

import map.ObjectSet;
import parse.Parse;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaReferenceGeometry.Hit;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Tools;
import structures.ByteBuilder;

/** Compares full legacy-GA and captured-batch evidence without repeating searches.
 * Both saved generations are validated against the same catalog and model roster.
 * Old native tables represented by saved hits are read once each: physical rows,
 * coordinates, printed scores and --cut_ga footers must match before their !/?
 * markers are used. Old zero-hit inputs remain covered by completed saved metadata.
 * Legacy saved GA_pass is never used as native inclusion. Historical model-file
 * provenance remains the accepted generation's responsibility, not a new file hash.
 * Sorted coordinate multisets retain input/target/model/strand/sequence bounds/
 * model bounds/truncation. Report paths and lines are witnesses, never join keys.
 * Raw retention, native-included retention and equal printed-score multiplicity
 * are independent comparisons. No paired unrounded-score change is inferred.
 * State is confined to one reporting thread; retained event objects come from
 * the shared loader, with primitive flag/counter arrays and reference-array sorts.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceRetention{

	static void process(final String input, final String previous, final String output,
		final NcrnaEvidenceCatalog catalog, final Roster roster){
		final NcrnaReferenceGeometry old=NcrnaReferenceGeometry.load(previous, catalog, roster, true);
		final NcrnaReferenceGeometry current=NcrnaReferenceGeometry.load(input, catalog, roster, true);
		require(!old.profiles().batchEvidence() && current.profiles().batchEvidence(),
			"Retention requires legacy whole-input controls and original-GA-bound captured batch candidates");
		require(old.profiles().observed()==catalog.size() && current.profiles().observed()==catalog.size(),
			"Both retention generations must cover the same complete catalog");
		final Hit[] before=old.hits().clone(), after=current.hits().clone();
		final boolean[] included=readNative(before, roster);
		Arrays.sort(before, SCORE_ORDER);
		Arrays.sort(after, SCORE_ORDER);
		compare(before, after, included, roster, output);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Original Native Rows    ----------------*/
	/*--------------------------------------------------------------*/

	/** Re-reads actual flags and verifies every saved event against its physical native row. */
	private static boolean[] readNative(final Hit[] hits, final Roster roster){
		Arrays.sort(hits, NATIVE_ORDER);
		final boolean[] included=new boolean[hits.length];
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		final RrnaCmTruthBuilder.Fields f=new RrnaCmTruthBuilder.Fields();
		final ByteBuilder target=new ByteBuilder();
		for(int first=0; first<hits.length; ){
			int end=first+1;
			while(end<hits.length && hits[end].input==hits[first].input && Tools.equals(hits[end].table, hits[first].table)){end++;}
			final String path=decode(hits[first].table);
			require(new File(path).isAbsolute() && seen.add(path), "Each legacy input must retain a distinct absolute native report");
			int next=first, programs=0, versions=0, commands=0, completed=0;
			long line=0;
			try(Input in=new Input(path)){
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					line++;
					if(row.length==0){continue;}
					if(completed!=0){throw new IllegalArgumentException("Content follows native control completion: "+path);}
					if(row[0]=='#'){
						final String text=new String(row, StandardCharsets.UTF_8);
						if(text.startsWith("# Program:")){
							require(text.substring(10).trim().equals("cmsearch"), "Expected native cmsearch control");
							programs++;
						}else if(text.startsWith("# Version:")){
							require(text.substring(10).trim().startsWith("1.1.5 "), "Expected captured Infernal1.1.5 control format");
							versions++;
						}else if(text.startsWith("# Option settings:")){
							checkCommand(text.substring(18).trim());
							commands++;
						}else if(text.equals("# [ok]")){completed++;}
						continue;
					}
					if(next>=end || hits[next].line!=line){
						throw new IllegalArgumentException("Native control has a missing/extra/reordered physical event: "+path+':'+line);
					}
					final Hit hit=hits[next++];
					final Model model=roster.query(hit.query);
					f.set(row);
					require(f.eq(4, "cm") && f.eq(3, model.accession) && f.eq(2, model.queryName()), "Native control model differs from its saved event");
					BuildNcrnaReference.field(target.clear(), row, f.start[0], f.end[0]);
					require(same(target, hit.target.bytes), "Native control target differs from its saved literal target");
					final int from=f.number(7), to=f.number(8);
					require(Math.min(from, to)==hit.start && Math.max(from, to)==hit.stop
						&& (hit.strand=='+' ? from<=to && f.eq(9, "+") : from>=to && f.eq(9, "-"))
						&& f.number(5)==hit.from && f.number(6)==hit.to && f.eq(10, TRUNCATIONS[hit.truncation]),
						"Native control alignment differs from its saved coordinate/truncation identity");
					final double bits=Parse.parseDouble(row, f.start[14], f.end[14]);
					require(Double.isFinite(bits) && bits==hit.bits, "Native control printed score differs from saved evidence");
					require(f.eq(16, "!") || f.eq(16, "?"), "Unknown native control inclusion marker");
					included[hit.index]=f.eq(16, "!");
				}
			}
			require(next==end && programs==1 && versions==1 && commands==1 && completed==1,
				"Native control rows/program/version/original-GA command/completion must all be present: "+path);
			first=end;
		}
		return included;
	}

	/** Only original GA controls may supply the previous generation's native flags. */
	private static void checkCommand(final String command){
		int ga=0;
		for(final String token:command.split(" +")){
			if(token.equals("--cut_ga")){ga++;}
			require(!token.equals("--cut_nc") && !token.equals("--cut_tc") && !token.equals("--incT")
				&& !token.equals("--incE") && !token.equals("-T") && !token.equals("-E"),
				"Original-GA control cannot override reporting or inclusion thresholds");
		}
		require(ga==1, "Original-GA control requires exactly one --cut_ga option");
	}

	/** Inverts only the four escapes emitted by BuildNcrnaReference.field. */
	private static String decode(final byte[] encoded){
		final ByteBuilder b=new ByteBuilder(encoded.length);
		for(int i=0; i<encoded.length; i++){
			if(encoded[i]!='%'){b.append(encoded[i]); continue;}
			require(i+2<encoded.length, "Incomplete saved path escape");
			final int a=Character.digit((char)encoded[++i], 16), z=Character.digit((char)encoded[++i], 16);
			require(a>=0 && z>=0, "Invalid saved path escape");
			final int value=a*16+z;
			require(value=='%' || value=='\t' || value=='\r' || value=='\n', "Unrecognized saved path framing escape");
			b.append((byte)value);
		}
		return new String(b.array, 0, b.length(), StandardCharsets.UTF_8);
	}
	private static boolean same(final ByteBuilder b, final byte[] expected){
		if(b.length()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(b.array[i]!=expected[i]){return false;}}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Multiset Comparison   ----------------*/
	/*--------------------------------------------------------------*/

	private static void compare(final Hit[] before, final Hit[] after, final boolean[] oldIncluded,
		final Roster roster, final String output){
		final long[][] totals=new long[roster.queried()+1][COLUMNS];
		final long[] row=new long[COLUMNS];
		final ByteBuilder b=new ByteBuilder();
		int i=0, j=0;
		try(Output out=new Output(output+SUFFIXES[0])){
			out.writer.println("input_id\ttarget\taccession\tstrand\tstart\tstop\tmodel_from\tmodel_to\tnative_truncation\t"+COUNTS
				+"\told_example_table\told_example_line\tnew_example_table\tnew_example_line");
			while(i<before.length || j<after.length){
				final Hit key=j==after.length || i<before.length && coordinate(before[i], after[j])<=0 ? before[i] : after[j];
				int ie=i, je=j;
				while(ie<before.length && coordinate(key, before[ie])==0){ie++;}
				while(je<after.length && coordinate(key, after[je])==0){je++;}
				Arrays.fill(row, 0);
				row[0]=ie-i;
				row[1]=je-j;
				row[2]=Math.min(row[0], row[1]);
				row[3]=row[0]-row[2];
				row[4]=row[1]-row[2];
				for(int k=i; k<ie; k++){if(oldIncluded[before[k].index]){row[5]++;}}
				for(int k=j; k<je; k++){if(after[k].included){row[6]++;}}
				row[7]=Math.min(row[5], row[6]);
				row[8]=row[5]-row[7];
				row[9]=row[6]-row[7];
				final long oldExcluded=row[0]-row[5], newExcluded=row[1]-row[6];
				final long excludedRetained=Math.min(oldExcluded, newExcluded);
				row[14]=Math.min(row[8], newExcluded-excludedRetained);
				row[15]=Math.min(oldExcluded-excludedRetained, row[9]);
				row[10]=row[9]-row[15];
				row[11]=newExcluded-excludedRetained-row[14];
				int a=i, z=j;
				while(a<ie && z<je){
					final int order=Double.compare(before[a].bits, after[z].bits);
					if(order==0){
						row[12]++;
						a++;
						z++;
					}
					else if(order<0){a++;}else{z++;}
				}
				row[13]=row[2]-row[12];
				conserve(row);
				for(int k=0; k<COLUMNS; k++){
					totals[key.query][k]+=row[k];
					totals[roster.queried()][k]+=row[k];
				}
				if(row[3]+row[4]+row[8]+row[9]+row[13]>0){
					b.clear().append(key.input.genome.id).tab().append(key.target.bytes).tab().append(roster.query(key.query).accession).tab();
					b.append(key.strand).tab().append(key.start).tab().append(key.stop).tab().append(key.from).tab().append(key.to).tab().append(key.truncation);
					for(final long value:row){b.tab().append(value);}
					example(b, i<ie ? before[i] : null);
					example(b, j<je ? after[j] : null);
					out.writer.print(b.nl());
				}
				i=ie;
				j=je;
			}
		}
		require(totals[roster.queried()][0]==before.length && totals[roster.queried()][1]==after.length,
			"Full retention summary must conserve both saved event populations");
		try(Output out=new Output(output+SUFFIXES[1])){
			out.writer.println("accession\t"+COUNTS+"\tstatus");
			for(int model=0; model<totals.length; model++){
				conserve(totals[model]);
				b.clear().append(model==roster.queried() ? "ALL_QUERIED" : roster.query(model).accession);
				for(final long value:totals[model]){b.tab().append(value);}
				out.writer.print(b.append("\tRAW_NATIVE_RETENTION_NOT_ADJUDICATED\n"));
			}
		}
	}

	/** Same-flag matches take priority when classifying raw additions in duplicate groups. */
	private static void conserve(final long[] r){
		for(final long value:r){require(value>=0, "Multiset counts cannot be negative");}
		require(r[0]==r[2]+r[3] && r[1]==r[2]+r[4] && r[5]==r[7]+r[8] && r[6]==r[7]+r[9]
			&& r[4]==r[10]+r[11] && r[2]==r[12]+r[13], "Raw/native-included/printed-score multiset accounting must conserve every event");
	}
	/** A representative member of the coordinate group, not a uniquely assigned lost duplicate. */
	private static void example(final ByteBuilder b, final Hit hit){
		if(hit==null){b.append("\tNA\t-1");}
		else{b.tab().append(hit.table).tab().append(hit.line);}
	}
	private static int coordinate(final Hit a, final Hit b){
		int c=Integer.compare(a.input.index, b.input.index);
		if(c==0){c=Tools.compare(a.target.bytes, b.target.bytes);}
		if(c==0){c=Integer.compare(a.query, b.query);}
		if(c==0){c=Character.compare(a.strand, b.strand);}
		if(c==0){c=Integer.compare(a.start, b.start);}
		if(c==0){c=Integer.compare(a.stop, b.stop);}
		if(c==0){c=Integer.compare(a.from, b.from);}
		if(c==0){c=Integer.compare(a.to, b.to);}
		return c==0 ? Integer.compare(a.truncation, b.truncation) : c;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private static final Comparator<Hit> SCORE_ORDER=(a, b)->{
		final int c=coordinate(a, b);
		return c==0 ? Double.compare(a.bits, b.bits) : c;
	};
	private static final Comparator<Hit> NATIVE_ORDER=(a, b)->{
		int c=Integer.compare(a.input.index, b.input.index);
		if(c==0){c=Tools.compare(a.table, b.table);}
		return c==0 ? Long.compare(a.line, b.line) : c;
	};
	private static final String[] TRUNCATIONS={"no", "5'", "3'", "5'&3'"};
	private static final int COLUMNS=16;
	private static final String COUNTS="old_raw\tnew_raw\traw_retained\traw_lost\traw_added\told_native_included\tnew_native_included\tincluded_retained\tincluded_lost\tincluded_gained\traw_added_included\traw_added_excluded\tequal_printed_score_matches\tcoordinate_matches_with_changed_printed_score\tincluded_to_excluded\texcluded_to_included";
	static final String[] SUFFIXES={".retention_differences.tsv.gz", ".retention_summary.tsv"};
}
