package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import fileIO.ReadWrite;
import prot.DigestSuffix;
import prot.FamilyAlignmentMetric;
import prot.FamilyCalibration;
import prot.FamilyHeuristicMetadata;
import prot.FamilyGateCodec;
import prot.FamilyMemberPreparation;
import prot.FamilyMemberPreparation.Prepared;
import prot.FamilyNegativeUpperCalibration;
import prot.FamilyNegativeUpperCalibration.MetricResult;
import prot.FamilyNegativeUpperCalibration.MetricWorker;
import prot.FamilyNegativeUpperCalibration.NegativePanel;
import prot.FamilyNegativeUpperCalibration.PackagedCalibration;
import prot.FamilyTailCascade;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyReferenceExport;
import prot.FamilyReferenceLoader;
import prot.MemberLengthRange;
import shared.Shared;
import structures.ByteBuilder;

/** Literal selection, denominator and absent-evidence checks for the native negative worker.
 * @author Brian Bushnell, Keqing
 */
public final class FamilyNegativeMetricTest{

	public static void main(String[] args)throws Exception{
		check(args.length==1,"Supply a fresh fixture directory");
		final Path root=Paths.get(args[0]); check(!Files.exists(root),"Do not overwrite negative fixture evidence"); Files.createDirectory(root);
		Shared.AMINO_IN=true; Shared.setThreads(1); ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final String core=repeat('A',40), shortCore=repeat('A',20);
		final Path negatives=root.resolve("negatives.fa"), medians=root.resolve("medians.tsv"), hydro=root.resolve("hydro.tsv");
		write(negatives,">q0\n"+core+"\n>q1\nR"+core+"\n>q2\n"+shortCore+"\n>q3\n"+repeat('R',40)+"\n");
		write(medians,"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n"+
			"normal\t1\t40\t40\t40\t1\nshortref\t1\t40\t40\t40\t1\n");
		write(hydro,"family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tnormal\t1\t1.8\n1\tshortref\t1\t1.8\n");
		final FamilyAlignmentMetric metric=new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,28);
		final FamilyMemberPreparation prep=preparer(metric);
		final Prepared normal=prep.prepare(0,"normal",new byte[][]{bytes(core)},bytes(core));
		final Prepared shortref=prep.prepare(1,"shortref",new byte[][]{bytes(core)},bytes(shortCore));
		final FamilyHeuristicMetadata metadata=new FamilyHeuristicMetadata(new Prepared[]{normal,shortref},medians,
			DigestSuffix.file(medians.toString()),hydro,DigestSuffix.file(hydro.toString()));
		final FamilyTailCascade.ScalarBounds bounds=new FamilyTailCascade.ScalarBounds(
			new double[]{-4.5,-1,0,0,0},new double[]{4.5,1,1,1,20},1);
		final NegativePanel panel=FamilyNegativeUpperCalibration.loadPreparationPanel(negatives,4,20,41,bounds,
			ProteinKmerGeometry.exact(5),ProteinPositionMode.PLAIN,1,0,null);
		final PackagedCalibration population=new PackagedCalibration(new Prepared[]{normal,shortref},panel,metric,metadata,
			new int[]{0,1},new int[]{1,1});
		final MetricWorker worker=population.newMetricWorker(4);
		final int[] observed={0};
		final MetricResult result=worker.scoreFamily(0,(f,q,h,c,a)->{
			check(f==0 && q!=3 && h>=1,"Raw-hit admission precedes selected-pair observation"); observed[0]++;
			if(q==1){check(a.measured() && a.identityPct()==100 && a.blosum()==1,"Leading R is clipped with aligned-query normalization");}
		});
		check(result.scanned==4 && result.admitted==3 && result.selected==3 && observed[0]==3,
			"Three admitted queries are selected from four scanned; top is a limit, not invented observations");
		check(result.eligible==2 && result.ineligible()==1 && result.count(FamilyAlignmentMetric.State.INPUT_TOO_SHORT)==1,
			"The twenty-residue query is counted but supplies no normalized upper");
		check(result.identityUpper()==100 && result.blosumUpper()==1 && result.metric.same(metric),"Detached maxima retain the measured metric");
		check(result.familyRank==0 && result.familyId.equals("normal"),"Export evidence keeps physical family identity beyond batch order");
		check(panel.rawCosineMaximum(normal)==1,"All-panel raw cosine remains independently measured");
		final MetricResult absent=worker.scoreFamily(1,null);
		check(absent.selected==3 && absent.eligible==0 && absent.count(FamilyAlignmentMetric.State.INPUT_TOO_SHORT)==3,
			"A short reference makes every selected alignment unavailable");
		reject(()->absent.identityUpper(),"Absent identity upper must throw");
		reject(()->absent.blosumUpper(),"Absent BLOSUM upper must throw");
		check(result.eligible==2 && result.count(FamilyAlignmentMetric.State.INPUT_TOO_SHORT)==1,"Later worker calls cannot mutate earlier result counts");
		final MetricResult none=new PackagedCalibration(new Prepared[]{normal},panel,metric,metadata,new int[]{0},new int[]{1024})
			.newMetricWorker(4).scoreFamily(0,null);
		check(none.scanned==4 && none.admitted==0 && none.selected==0 && !none.measured(),"Zero admission is explicit absent evidence");
		reject(()->none.blosumUpper(),"Zero selected is not empirical upper zero");
		final int[] winner={-1};
		population.newMetricWorker(1).scoreFamily(0,(f,q,h,c,a)->winner[0]=q);
		check(winner[0]==0,"Exact forty-A query outranks length-mismatched prefixes under the released pivot");
		reject(()->population.newWorker(4),"Legacy result consumers must reject explicit alignment states");
		reject(()->new PackagedCalibration(new Prepared[]{normal},panel),"Legacy construction must reject new member metric");
		reject(()->new PackagedCalibration(new Prepared[]{normal},panel,metric,metadata,new int[]{1},new int[]{1}),"Batch map must match family identity and rank");
		reject(()->worker.scoreFamily(0,(f,q,h,c,a)->{throw new IllegalArgumentException("observer abort");}),"Observer failures propagate");
		check(worker.scoreFamily(0,null).eligible==2,"Worker resets after observer failure");
		final Prepared legacy=preparer(FamilyAlignmentMetric.LEGACY).prepare(0,"normal",new byte[][]{bytes(core)},bytes(core));
		final FamilyNegativeUpperCalibration.FamilyResult old=new PackagedCalibration(new Prepared[]{legacy},panel).newWorker(4).scoreFamily(0,null);
		check(old.aligned==4 && old.observations==4 && old.identity==100 && old.blosum==1 && old.cosineMax==result.cosineMax,
			"Legacy forced-all selection and measured maxima remain intact");
		final FamilyReference measuredExport=FamilyReferenceExport.fromPrepared(normal,1,new MetricResult[]{none,result},0,0);
		check(measuredExport.alignmentCalibration==FamilyReference.AlignmentCalibration.ENABLED &&
			measuredExport.identityUpper==Math.nextUp(110f) && measuredExport.blosumFinal==FamilyGateCodec.roundUpper(1.1),
			"Only measured panel maxima enter the original margin formula");
		final FamilyReference explicitDefault=FamilyReferenceExport.fromPrepared(normal,1,new MetricResult[]{none,result},0,0,
			FamilyGateCodec.DEFAULT_IDENTITY_UPPER_MULTIPLIER,FamilyGateCodec.DEFAULT_IDENTITY_UPPER_ADDITIVE);
		check(explicitDefault.identityUpper==measuredExport.identityUpper && explicitDefault.blosumFinal==measuredExport.blosumFinal,
			"Explicit default padding must preserve default metric export gates");
		check(FamilyReferenceExport.fromPrepared(normal,1,new MetricResult[]{result},0,0,1,1).identityUpper==101,
			"Metric export can use +1 identity padding from raw measured maxima");
		final FamilyReference noUpper=FamilyReferenceExport.fromPrepared(normal,1,new MetricResult[]{none},0,0);
		check(noUpper.identityUpper==100 && noUpper.blosumFinal==1,"Absent negative evidence omits upper margin terms and retains measured member bounds");
		check(FamilyReferenceExport.fromPrepared(normal,1,new MetricResult[]{none},1.2,0).blosumFinal==FamilyGateCodec.roundUpper(1.2),
			"A configured global floor remains active when negative evidence is unavailable");
		final FamilyReference disabled=FamilyReferenceExport.fromPrepared(shortref,1,new MetricResult[]{absent},0,1.5);
		check(disabled.alignmentCalibration==FamilyReference.AlignmentCalibration.UNAVAILABLE && disabled.identityLower==0 &&
			disabled.identityUpper==0 && disabled.blosumLower==0 && disabled.blosumFinal==0 && disabled.cosineAccept>1,
			"Absent member alignments disable only alignment channels, retaining real cosine gates");
		reject(()->FamilyReferenceExport.fromPrepared(normal,1,new MetricResult[]{none},0,0,Double.NaN,0),
			"Metric export must reject invalid padding even without measured upper evidence");
		reject(()->FamilyReferenceExport.fromPrepared(shortref,1,new MetricResult[]{absent},0,0,0.5,-1),
			"Metric export must reject invalid padding even when alignment gates are unavailable");
		final FamilyReference projected=FamilyReferenceExport.projectMonomers(disabled);
		check(projected.alignmentCalibration==FamilyReference.AlignmentCalibration.UNAVAILABLE && projected.spectrum.length==20,
			"Projection must preserve unavailable alignment state");
		reject(()->FamilyReferenceCodec.appendFiveLines(projected,new ByteBuilder()),"Projection/scalar copies preserve explicit framing and reject legacy serialization");
		final FamilyReference projectedMeasured=FamilyReferenceExport.projectMonomers(measuredExport);
		final FamilyReference projectedExplicitDefault=FamilyReferenceExport.projectMonomers(explicitDefault);
		final FamilyReference fallback=copyIdentity(projectedMeasured,"fallback",2,projectedMeasured.identityUpper);
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(FamilyReferenceCodec.ALIGNMENT_MONOMER_SCHEMA,
			5,1024,0,3,3,0,FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING,0,metric);
		final ByteBuilder explicitBytes=new ByteBuilder(),defaultBytes=new ByteBuilder();
		FamilyReferenceCodec.appendFiveLines(projectedMeasured,defaultBytes,header);
		FamilyReferenceCodec.appendFiveLines(projectedExplicitDefault,explicitBytes,header);
		check(defaultBytes.toString().equals(explicitBytes.toString()),
			"Explicit default padding must serialize identically to implicit default gates and payloads");
		final ByteBuilder references=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(header,references);
		FamilyReferenceCodec.appendFiveLines(projectedMeasured,references,header);
		FamilyReferenceCodec.appendFiveLines(projected,references,header);
		FamilyReferenceCodec.appendFiveLines(fallback,references,header);
		final Path parent=root.resolve("references.tsv"),raw=root.resolve("raw_calibration.tsv"),
			rawMetricMismatch=root.resolve("raw_metric_mismatch.tsv"),rawRosterMismatch=root.resolve("raw_roster_mismatch.tsv"),
			padded=root.resolve("padded"),paddedFractional=root.resolve("padded_fractional");
		final Path parentAbs=parent.toAbsolutePath().normalize(),rawAbs=raw.toAbsolutePath().normalize(),
			paddedAbs=padded.toAbsolutePath().normalize();
		write(parent,references.toString());
		write(raw,rawCalibration(metric.name(),28,"fallback","NA"));
		write(rawMetricMismatch,rawCalibration(metric.name(),27,"fallback","NA"));
		write(rawRosterMismatch,rawCalibration(metric.name(),28,"wrong","NA"));
		rejectRun(new String[]{"rawcalibration="+rawAbs,"rawcalibrationsha80=unused"},
			"Raw-calibration replay inputs require explicit identitypadding mode");
		rejectRun(new String[]{"mode=identitypadding","parent="+parentAbs,"parentsha80=unused",
			"rawcalibration="+rawAbs,"out="+root.resolve("sha_rejected").toAbsolutePath().normalize()},
			"Identity-padding replay must reject routine input pins");
		rejectRun(new String[]{"mode=identitypadding","parent="+parentAbs,
			"rawcalibration="+rawMetricMismatch.toAbsolutePath().normalize(),
			"out="+root.resolve("metric_mismatch").toAbsolutePath().normalize()},
			"Saved raw metric columns must match the parent reference header");
		rejectRun(new String[]{"mode=identitypadding","parent="+parentAbs,
			"rawcalibration="+rawRosterMismatch.toAbsolutePath().normalize(),
			"out="+root.resolve("roster_mismatch").toAbsolutePath().normalize()},
			"Saved raw rows must match the parent reference roster exactly");
		runIdentityPadding(parentAbs,rawAbs,paddedAbs,1,1);
		runIdentityPadding(parentAbs,rawAbs,paddedFractional.toAbsolutePath().normalize(),1,0.25);
		final FamilyReference[] replayed=FamilyReferenceLoader.loadMetricBound(padded.resolve("references.tsv").toString(),1,10).all;
		final FamilyReference[] fractional=FamilyReferenceLoader.loadMetricBound(paddedFractional.resolve("references.tsv").toString(),1,10).all;
		check(replayed[0].identityUpper==101 && replayed[0].blosumFinal==projectedMeasured.blosumFinal,
			"Saved raw identity maxima can be replayed with a new identity-only margin");
		check(fractional[0].identityUpper==FamilyGateCodec.upper(100,100,1,0.25),
			"Saved raw identity maxima can be replayed with a fractional identity padding margin");
		check(replayed[1].alignmentCalibration==FamilyReference.AlignmentCalibration.UNAVAILABLE && replayed[1].identityUpper==0,
			"Saved unavailable alignment rows stay explicit and inert through identity padding replay");
		check(replayed[2].identityUpper==100,
			"Saved rows without negative identity evidence fall back to the raw lower and thirty-percent floor");
		check(!Files.exists(padded.resolve("inputs.tsv")) && !Files.exists(padded.resolve("result.sha80.tsv")),
			"Identity-padding replay validates structure instead of publishing routine input/output hashes");
		final ByteBuilder restored=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(header,restored);
		FamilyReferenceCodec.appendFiveLines(copyIdentity(replayed[0],projectedMeasured.id,projectedMeasured.rank,
			projectedMeasured.identityUpper),restored,header);
		FamilyReferenceCodec.appendFiveLines(replayed[1],restored,header);
		FamilyReferenceCodec.appendFiveLines(copyIdentity(replayed[2],fallback.id,fallback.rank,fallback.identityUpper),
			restored,header);
		check(references.toString().equals(restored.toString()),
			"Identity replay must preserve the header and every non-identity payload byte");
		reject(()->FamilyReferenceExport.fromPrepared(normal,1,new MetricResult[]{absent},0,0),"Upper evidence from another physical family cannot be exported");
		reject(()->FamilyReferenceExport.fromPrepared(normal,1,new MetricResult[0],0,0),"Missing panel results are not unavailable evidence");
		reject(()->FamilyReferenceExport.fromPrepared(legacy,1,new MetricResult[]{result},0,0),"Legacy members cannot consume new-metric uppers");
		reject(()->FamilyReferenceExport.applyAcceptOverrides(root.resolve("unread.tsv"),new FamilyReference[]{noUpper}),
			"Legacy override tables cannot silently change new-metric gates");
		System.out.println("FAMILY_NEGATIVE_METRIC_TEST_PASS admission=true clipping=true accounting=true unavailable=true legacy=true");
		System.out.println("FAMILY_METRIC_EXPORT_TEST_PASS evidence_bound=true unavailable=true copies=true legacy_guard=true");
	}
	private static FamilyMemberPreparation preparer(FamilyAlignmentMetric metric){
		return new FamilyMemberPreparation(ProteinKmerGeometry.exact(5),4,2,null,FamilyCalibration.MemberRetention.METRIC_SPECIFIC,
			false,new MemberLengthRange(40,99),ProteinPositionMode.FIXED,metric);
	}
	private static String rawCalibration(String metricName,int minimumSpan,String fallbackId,String fallbackUpper){
		return "roster_rank\tfamily_id\tmembers\tcosine_lower\tcosine_retained\tidentity_lower\tidentity_retained\tblosum_lower\tblosum_retained\tcosine_upper\tidentity_upper\tblosum_upper\talignment_mode\talignment_min_span\tmember_alignment_measured\tmember_alignment_failed\tidentity_quota_achievable\tblosum_quota_achievable\n"+
			"0\tnormal\t1\t1\t1\t100\t1\t1\t1\t1\t100\t1\t"+metricName+"\t"+minimumSpan+"\t1\t0\ttrue\ttrue\n"+
			"1\tshortref\t1\t1\t1\tNA\tNA\tNA\tNA\t1\t88\t0.5\t"+metricName+"\t"+minimumSpan+"\t0\t1\tfalse\tfalse\n"+
			"2\t"+fallbackId+"\t1\t1\t1\t100\t1\t1\t1\t1\t"+fallbackUpper+"\t1\t"+metricName+"\t"+minimumSpan+"\t1\t0\ttrue\ttrue\n";
	}
	private static void runIdentityPadding(Path parent,Path raw,Path out,double multiplier,double additive)throws Exception{
		FamilyReferenceExport.run(new String[]{"mode=identitypadding","parent="+parent,
			"rawcalibration="+raw,"identityuppermultiplier="+multiplier,
			"identityupperadditive="+additive,"out="+out,"threads=1"});
	}
	private static FamilyReference copyIdentity(FamilyReference f,String id,int rank,float identityUpper){
		return new FamilyReference(id,rank,f.members,f.residues,f.validWords,new prot.ProteinSequence(id,f.consensusRaw.clone()),
			f.ddl.clone(),f.spectrum.clone(),f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,identityUpper,
			f.blosumFinal,f.consensusRaw.clone(),f.scalarMask,f.scalarBounds.clone(),f.alignmentCalibration);
	}
	private static String repeat(char acid,int n){final char[] out=new char[n];Arrays.fill(out,acid);return new String(out);}
	private static byte[] bytes(String value){return value.getBytes(StandardCharsets.US_ASCII);}
	private static void write(Path path,String value)throws Exception{Files.write(path,bytes(value));}
	private static void check(boolean valid,String why){if(!valid){throw new AssertionError(why);}}
	private static void reject(Runnable action,String why){
		try{action.run();}catch(IllegalArgumentException|IllegalStateException expected){return;}
		throw new AssertionError(why);
	}
	private static void rejectRun(String[] args,String why)throws Exception{
		try{FamilyReferenceExport.run(args);}catch(IllegalArgumentException|IllegalStateException expected){return;}
		throw new AssertionError(why);
	}
}
