package test.prot;

import java.util.Arrays;
import java.util.Random;

import ddl.DDLFamilyProfiles;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import prot.FamilyAlignmentMetric;
import prot.FamilyAlignmentMetric.State;
import prot.FamilyCalibration;
import prot.FamilyMemberPreparation;
import prot.FamilyMemberPreparation.Prepared;
import prot.FamilyReferenceExport;
import prot.MemberLengthRange;
import structures.ByteBuilder;

/** Literal member geometry and full-population retention controls; no scientific input panels.
 * @author Brian Bushnell, Keqing
 */
public final class FamilyMemberMetricPreparationTest{

	public static void main(String[] args){
		check(args.length==0,"Member metric fixture takes no file arguments");
		final FamilyAlignmentMetric pinned=new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,28);
		final FamilyMemberPreparation worker=worker(pinned);
		final byte[] core=repeat('A',40),extended=repeat('A',41); extended[0]='R';
		final byte[][] members={core,extended};
		final Prepared legacy=worker(FamilyAlignmentMetric.LEGACY).prepare(0,"family",members,core);
		final Prepared clipped=worker.prepare(0,"family",members,core);
		check(legacy.lower!=null && legacy.measuredLower==null && clipped.lower==null && clipped.measuredLower!=null,
			"A new metric must not masquerade as legacy all-member quantiles");
		check(legacy.lower[2].cutoff==156.0/165 && legacy.lower[1].cutoff==100.0*40/41,
			"Historical full-query scores retain the literal leading-R penalty");
		check(clipped.measuredLower[1].cutoff()==100 && clipped.measuredLower[2].cutoff()==1 &&
			clipped.measuredLower[1].retained==2 && clipped.measuredLower[2].retained==2,
			"Clipped forty-A cores retain both members at identity100 and ratio1");
		check(clipped.alignmentMetric.same(pinned) && clipped.alignmentCounts[State.MEASURED.ordinal()]==2,
			"Prepared output preserves its metric and measured population");
		check(Arrays.equals(legacy.plain,clipped.plain) && Arrays.equals(legacy.positional,clipped.positional) &&
			Arrays.equals(legacy.spectrum.counts,clipped.spectrum.counts) && Arrays.equals(legacy.statistics,clipped.statistics) &&
			legacy.lower[0].cutoff==clipped.measuredLower[0].cutoff(),
			"Changing alignment leaves member sketches, spectra, scalars and cosine arithmetic unchanged");
		check(FamilyReferenceExport.fromPrepared(legacy,0,0,0,0,0).members==2,"Legacy prepared export still works");
		fails(()->FamilyReferenceExport.fromPrepared(clipped,0,0,0,0,0),
			"Unbound legacy negative uppers must not be exported with new-metric member gates");

		final byte[] negative=repeat('R',40);
		for(int n:new int[]{100,400}){
			final byte[][] population=new byte[n][]; Arrays.fill(population,core); population[n-1]=negative;
			final Prepared mixed=worker.prepare(1,"mixed",population,core);
			check(mixed.alignmentCounts[State.MEASURED.ordinal()]==n-1 && mixed.alignmentCounts[State.NO_ALIGNMENT.ordinal()]==1,
				"Every assigned member has exactly one measured or failed alignment outcome");
			check(mixed.measuredLower[0].measured==n && mixed.measuredLower[0].failed==0,
				"A failed alignment does not remove its full-sequence cosine observation");
			for(int metric=1;metric<3;metric++){
				final FamilyCalibration.MemberCutoff cutoff=mixed.measuredLower[metric];
				check(cutoff.assigned==n && cutoff.measured==n-1 && cutoff.failed==1 && cutoff.retained==n-1 &&
					cutoff.quotaAchievable==(n==400),"Unavailable members consume the original member-retention budget");
			}
		}
		final Prepared absent=worker.prepare(2,"absent",new byte[][]{negative},core);
		check(absent.measuredLower[0].measured() && !absent.measuredLower[1].measured() && !absent.measuredLower[2].measured(),
			"All-failed alignment retains a real cosine and explicit unavailable alignment cutoffs");
		fails(()->absent.measuredLower[1].cutoff(),"Unavailable identity cutoff cannot be read as zero");
		fails(()->absent.measuredLower[2].cutoff(),"Unavailable BLOSUM cutoff cannot be read as zero");
		final Prepared shortReference=worker.prepare(3,"short_ref",members,repeat('A',20));
		check(shortReference.alignmentCounts[State.INPUT_TOO_SHORT.ordinal()]==2 && !shortReference.measuredLower[1].measured(),
			"Reference span eligibility applies to every member without inventing alignment measurements");
		final byte[] shortCore=repeat('R',40),other=repeat('C',40); shortCore[39]=other[39]='A';
		final Prepared spanFailure=worker.prepare(4,"short_span",new byte[][]{shortCore},other);
		check(spanFailure.alignmentCounts[State.SPAN_TOO_SHORT.ordinal()]==1 && !spanFailure.measuredLower[2].measured(),
			"A positive one-A suffix is distinct from a qualified twenty-eight-residue alignment");
		final Prepared empty=worker.prepare(5,"empty",new byte[0][],core);
		check(empty.members==0 && empty.lower==null && empty.measuredLower==null && empty.alignmentCounts==null &&
			empty.alignmentMetric.same(pinned),"Zero-member output keeps the metric but no measured payload");
		check(worker.prepare(0,"family",members,core).measuredLower[2].cutoff()==1,
			"Worker reuse after failures and emptiness cannot retain previous score state");
		widths(pinned);
		System.out.println("FAMILY_MEMBER_METRIC_PREPARATION_TEST_PASS legacy=true clipping=true accounting=true guards=true");
	}

	/** Fold each original member independently, then count modal votes without using Profile. */
	private static void widths(FamilyAlignmentMetric metric){
		final byte[] alphabet="ACDEFGHIKLMNPQRSTVWY".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final Random random=new Random(20261106);
		final byte[][] members=new byte[5][80];
		for(int i=0;i<80;i++){members[0][i]=alphabet[random.nextInt(alphabet.length)];}
		for(int m=1;m<members.length;m++){
			members[m]=members[0].clone();
			for(int i=m;i<80;i+=13){members[m][i]=alphabet[random.nextInt(alphabet.length)];}
		}
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.exact(5);
		for(ProteinPositionMode mode:ProteinPositionMode.values()){
			final int partitions=mode.positional() ? 4 : 1,overlap=mode.positional() ? 2 : 0;
			final FamilyMemberPreparation old=new FamilyMemberPreparation(geometry,partitions,overlap,null,
				FamilyCalibration.MemberRetention.METRIC_SPECIFIC,true,new MemberLengthRange(40,99),mode,metric);
			final Prepared original=old.prepare(0,"width",members,members[0]);
			for(int width:new int[]{1024,512}){
				final FamilyMemberPreparation worker=new FamilyMemberPreparation(geometry,partitions,overlap,null,
					FamilyCalibration.MemberRetention.METRIC_SPECIFIC,true,new MemberLengthRange(40,99),mode,metric,width);
				final Prepared result=worker.prepare(0,"width",members,members[0]);
				check(result.width==width && result.plain.length==width && result.positional.length==width &&
					result.positionalSupport.length==width && result.positionalOccupancy.length==width,
					"Both modal rows and their evidence arrays retain the explicitly requested width");
				check(Arrays.equals(original.statistics,result.statistics) && Arrays.equals(original.spectrum.counts,result.spectrum.counts) &&
					Arrays.equals(original.alignmentCounts,result.alignmentCounts) && original.validWords==result.validWords,
					"Sketch width cannot change member scalars, spectra, window counts or alignment outcomes");
				for(int metricIndex=0;metricIndex<3;metricIndex++){
					check(original.measuredLower[metricIndex].cutoff()==result.measuredLower[metricIndex].cutoff(),
						"Member-derived alignment and cosine gates are independent of sketch width");
				}
				if(width==1024){
					check(Arrays.equals(original.plain,result.plain) && Arrays.equals(original.positional,result.positional) &&
						Arrays.equals(original.positionalSupport,result.positionalSupport) &&
						Arrays.equals(original.positionalOccupancy,result.positionalOccupancy),
						"Existing constructor and explicit1024 preserve identical profile payloads");
				}else{
					foldedOracle(members,ProteinPositionMode.PLAIN,result.plain,null,null);
					foldedOracle(members,mode,result.positional,result.positionalSupport,result.positionalOccupancy);
				}
				final ByteBuilder row=new ByteBuilder();
				check(result.appendPositionalProfile(row),"A nonempty family emits its profile row");
				final String[] fields=row.toString().split("\t",-1);
				check(fields.length==10 && Integer.toString(width).equals(fields[2]),"Profile row declares its actual register width");
				for(int field=7;field<10;field++){
					check(fields[field].split(",",-1).length==width,"Serialized values/support/occupancy agree with the width column");
				}
				final Prepared empty=worker.prepare(1,"empty",new byte[0][],members[0]);
				check(empty.width==width && empty.plain==null && empty.positional==null &&
					!empty.appendPositionalProfile(row.clear()),"Empty families retain width metadata without invented votes");
				check(Arrays.equals(result.positional,worker.prepare(0,"width",members,members[0]).positional),
					"Worker reuse after emptiness must not retain another family's votes");
			}
		}
		for(int invalid:new int[]{0,256,513,2048}){
			fails(()->new FamilyMemberPreparation(geometry,4,2,null,FamilyCalibration.MemberRetention.METRIC_SPECIFIC,
				true,new MemberLengthRange(40,99),ProteinPositionMode.FIXED,metric,invalid),
				"Unsupported preparation widths fail before allocating member scratch");
		}
	}

	/** Exhaustive per-slot voting supplies an independent small-panel oracle, including unsigned-key ties. */
	private static void foldedOracle(byte[][] members,ProteinPositionMode mode,char[] actual,int[] support,int[] occupancy){
		final char[][] folded=new char[members.length][512];
		for(int m=0;m<members.length;m++){
			final char[] full=new char[1024];
			if(mode.positional()){
				DDLFamilyProfiles.sketchPositional(members[m],ProteinKmerGeometry.exact(5),0,4,2,mode,full,null);
			}else{DDLFamilyProfiles.sketch(members[m],5,0,full);}
			for(int b=0;b<512;b++){folded[m][b]=(char)Math.max(full[b],full[b+512]);}
		}
		for(int b=0;b<512;b++){
			char winner=0;int votes=0,occupied=0;
			for(char[] row:folded){
				final char candidate=row[b];if(candidate==0){continue;}occupied++;
				int count=0;for(char[] other:folded){if(other[b]==candidate){count++;}}
				if(count>votes || count==votes && candidate<winner){winner=candidate;votes=count;}
			}
			check(actual[b]==winner,"512 modal keys must be voted after per-member folding");
			if(support!=null){check(support[b]==votes && occupancy[b]==occupied,"Modal support and occupancy count original members");}
		}
	}

	private static FamilyMemberPreparation worker(FamilyAlignmentMetric metric){
		return new FamilyMemberPreparation(ProteinKmerGeometry.exact(5),4,2,null,
			FamilyCalibration.MemberRetention.METRIC_SPECIFIC,false,new MemberLengthRange(40,99),ProteinPositionMode.FIXED,metric);
	}
	private static byte[] repeat(char acid,int count){final byte[] raw=new byte[count]; Arrays.fill(raw,(byte)acid); return raw;}
	private static void fails(Runnable action,String why){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError(why);
	}
	private static void check(boolean valid,String why){if(!valid){throw new AssertionError(why);}}
}
