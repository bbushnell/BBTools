package test.prot;

import java.nio.charset.StandardCharsets;
import java.util.Random;

import parse.LineParser1;
import prot.FamilyGateCodec;
import structures.ByteBuilder;

/** Raw-bit and directional-rounding oracles for the shared A48 gate codec.
 * @author Sayu
 */
public final class FamilyGateCodecTest{

	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No codec test arguments");}
		for(int bits : new int[]{0,0x80000000,1,0x80000001,0x007fffff,0x00800000,0x3f800000,0x7f7fffff,0xff7fffff}){
			roundTrip(bits);
		}
		final Random random=new Random(4819321);
		for(int i=0; i<10000; i++){
			final int bits=random.nextInt(); final float value=Float.intBitsToFloat(bits);
			if(Float.isFinite(value)){roundTrip(bits); directed((double)value);}
		}
		for(double value : new double[]{0.3,-0.3,1.0/3,-1.0/3,8.896797153024911E-4,
			Double.MIN_VALUE,-Double.MIN_VALUE,(double)Float.MIN_VALUE/2,-(double)Float.MIN_VALUE/2,
			Math.nextDown(1.0),Math.nextUp(1.0),0.0,-0.0,Float.MAX_VALUE,-Float.MAX_VALUE}){directed(value);}
		for(double[] pair:new double[][]{{0,0},{30,30},{50,30},{30,36.79},{80,50},{100,100},{40,Math.nextDown(40/1.1)}}){
			final double bound=Math.max(pair[0],Math.max(30,Math.max(1.1*pair[1],pair[1]+2)));
			final float value=FamilyGateCodec.upper(pair[0],pair[1]);
			check(value>=bound && Math.nextDown(value)<bound,"Identity floor and margins use original doubles with one upward rounding");
		}
		check(FamilyGateCodec.upper(20,20)==30f && FamilyGateCodec.upper(50,30)==50f,"Floor and lower independently dominate identity acceptance");
		check(FamilyGateCodec.upper(30,36.79,1,1)==FamilyGateCodec.roundUpper(37.79),"Explicit identity additive margin is configurable");
		check(FamilyGateCodec.upper(30,50,1,1)==51f,"Unit identity multiplier disables proportional padding");
		reject(()->FamilyGateCodec.upper(0.3,Double.NaN));
		reject(()->FamilyGateCodec.upper(30,30,Double.NaN,1)); reject(()->FamilyGateCodec.upper(30,30,0.5,1));
		reject(()->FamilyGateCodec.upper(30,30,1,-1));
		reject(()->FamilyGateCodec.upper(-1,30)); reject(()->FamilyGateCodec.upper(30,101));
		final float sentinel=Float.intBitsToFloat(FamilyGateCodec.DISABLED_BITS);
		final ByteBuilder encoded=FamilyGateCodec.appendUpper(new ByteBuilder(),sentinel);
		check(FamilyGateCodec.disabled(FamilyGateCodec.readUpper(parser(encoded.toBytes()),0)),"Canonical sentinel bit roundtrip");
		check(sentinel==100f && Float.isFinite(sentinel),"Disabled identity upper is the finite percent endpoint");
		for(float percent:new float[]{0,50,Math.nextDown(100f),100f}){
			check(!(percent>sentinel),"No valid identity percent can pass a disabled shortcut");
		}
		check(FamilyGateCodec.readLower(parser(encoded.toBytes()),0)==100f,"Finite endpoint remains an ordinary finite value outside the identity-upper context");
		check(FamilyGateCodec.SCHEMA.equals("family_gate_a48_f32_v6"),"New margin policy requires a new schema");
		check(FamilyGateCodec.UPPER_HEADER.contains("\tidentity_pct_accept_a48f32_v3\t"),"Identity floor/margins have a versioned header");
		check(FamilyGateCodec.UPPER_HEADER.contains("\tcosine_accept_a48f32_v2\t"),"Final cosine must have a distinct header");
		for(double[] pair:new double[][]{{0,0},{0.3,0.3},{0.8,0.2},{0.1,0.8},{0.5,0.96},{0.5,0.9090909},{1,0},{1,1}}){
			final float value=FamilyGateCodec.cosineAccept(pair[0],pair[1]);
			final double raw=Math.max(pair[0],Math.max(0.72,Math.max(pair[1]*1.1,pair[1]+0.05)));
			check(Float.isFinite(value) && value>=raw && (double)Math.nextDown(value)<raw,"Cosine formula must round original doubles upward exactly once for every finite pair");
			final byte[] token=FamilyGateCodec.appendCosineAccept(new ByteBuilder(),value).toBytes();
			check(Float.floatToRawIntBits(FamilyGateCodec.readCosineAccept(parser(token),0))==Float.floatToRawIntBits(value),"Finite final cosine raw bits roundtrip");
		}
		check(FamilyGateCodec.cosineAccept(0.5,0.9090909)==1.0f,"Original-double cosine boundary is one, not a margin on a pre-rounded upper");
		for(double bad:new double[]{Double.NaN,Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY}){
			reject(()->FamilyGateCodec.cosineAccept(bad,0.5));
			reject(()->FamilyGateCodec.cosineAccept(0.5,bad));
		}
		check(FamilyGateCodec.UPPER_HEADER.endsWith("\tblosum_a_sq_final_a48f32_v4"),"Uncapped final formula must have a distinct header");
		for(double u : new double[]{Math.nextDown(-0.005),-0.005,Math.nextUp(-0.005),-0.001,-0.0,0.0,
			Math.nextDown(0.05),0.05,Math.nextUp(0.05)}){finalPair(-0.4,u);}
		check(FamilyGateCodec.finalBlosum(-0.4,-0.3)==0,"Negative upper floor");
		check(FamilyGateCodec.finalBlosum(-0.4,-0.001)>0,"Additive margin positive for upper in (-.005,0)");
		for(double[] pair : new double[][]{{0.3,0.3},{0.3,Math.nextDown(0.3)},{0.3,Math.nextUp(0.3)},
			{-0.3,-0.3},{-0.3,Math.nextUp(-0.3)},{0.0,-0.0},{Double.MIN_VALUE,0},
			{-Double.MIN_VALUE,-Float.MIN_VALUE},{Float.MAX_VALUE,Float.MAX_VALUE},{-Float.MAX_VALUE,-Float.MAX_VALUE},
			{0.1,0.11},{0.1,0.2},{0.2,0.1},{0.1,0.099},{-0.3,-0.29},
			{0.1,Math.nextDown(0.144/1.1)},{0.1,0.144/1.1},{0.1,Math.nextUp(0.144/1.1)},
			{Math.nextDown(0.144),0.2},{0.144,0.2},{Math.nextUp(0.144),0.2}}){
			finalPair(pair[0],pair[1]); finalPair(pair[1],pair[0]);
		}
		check(FamilyGateCodec.finalBlosum(0.1,0.11)>0.11,"Multiplier must raise an uncapped negative upper");
		check(FamilyGateCodec.finalBlosum(0.1,0.2)<FamilyGateCodec.finalBlosum(0.1,1),"Uncapped upper continues increasing");
		check(FamilyGateCodec.finalBlosum(0.1,0.2)>=0.22,"No0.144cap may erase a measured negative margin");
		check(FamilyGateCodec.finalBlosum(0.2,0.1)>=0.2,"Higher member lower must win");
		check(FamilyGateCodec.finalBlosum(0.1,0.099)>0.1,"Margin is applied before resolving a crossed raw pair");
		check((double)FamilyGateCodec.lower(0.3)<Math.nextDown(0.3),"Crossed fixture must expose downward-lower fallback risk");
		check((double)FamilyGateCodec.finalBlosum(0.3,Math.nextDown(0.3))>=0.3,"Crossed pair must not fall back to rounded-down lower");
		for(double bad : new double[]{Double.NaN,Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY,
			Math.nextUp((double)Float.MAX_VALUE),Math.nextDown(-(double)Float.MAX_VALUE)}){
			reject(()->FamilyGateCodec.lower(bad));
			reject(()->FamilyGateCodec.upper(bad,1));
			reject(()->FamilyGateCodec.finalBlosum(bad,1));
			reject(()->FamilyGateCodec.finalBlosum(0,bad));
			reject(()->FamilyGateCodec.upper(0,bad));
		}
		for(String bad : new String[]{"","00","01","0000000","/","p","-1","4ooooo","oooooo","1.0","Infinity"}){
			reject(()->FamilyGateCodec.readLower(parser(bad.getBytes(StandardCharsets.US_ASCII)),0));
			reject(()->FamilyGateCodec.readUpper(parser(bad.getBytes(StandardCharsets.US_ASCII)),0));
			reject(()->FamilyGateCodec.readFinalBlosum(parser(bad.getBytes(StandardCharsets.US_ASCII)),0));
			reject(()->FamilyGateCodec.readCosineAccept(parser(bad.getBytes(StandardCharsets.US_ASCII)),0));
		}
		for(int bits : new int[]{0x7f800000,0xff800000,0x7fc00000,0x7fc00001,0xffc00000,0x7f800001}){
			final float bad=Float.intBitsToFloat(bits);
			final byte[] token=new ByteBuilder().appendFloatA48(bad).toBytes();
			reject(()->FamilyGateCodec.readUpper(parser(token),0));
			reject(()->FamilyGateCodec.readLower(parser(token),0));
			reject(()->FamilyGateCodec.appendUpper(new ByteBuilder(),bad));
			reject(()->FamilyGateCodec.appendLower(new ByteBuilder(),bad));
			reject(()->FamilyGateCodec.readFinalBlosum(parser(token),0));
			reject(()->FamilyGateCodec.appendFinalBlosum(new ByteBuilder(),bad));
			reject(()->FamilyGateCodec.readCosineAccept(parser(token),0));
			reject(()->FamilyGateCodec.appendCosineAccept(new ByteBuilder(),bad));
		}
		final LineParser1 missing=parser("0".getBytes(StandardCharsets.US_ASCII));
		reject(()->FamilyGateCodec.readLower(missing,-1)); reject(()->FamilyGateCodec.readUpper(missing,1));
		reject(()->FamilyGateCodec.readFinalBlosum(missing,1));
		System.out.println("FAMILY_GATE_CODEC_TEST_PASS checks="+checks);
	}

	private static void roundTrip(int bits){
		final float value=Float.intBitsToFloat(bits);
		final ByteBuilder token=FamilyGateCodec.appendLower(new ByteBuilder(),value);
		check(token.length>=1 && token.length<=6,"Float A48 token length");
		final LineParser1 lp=parser(new ByteBuilder().append("prefix\t").append(token).append("\tsuffix").toBytes());
		check(Float.floatToRawIntBits(FamilyGateCodec.readLower(lp,1))==bits,"Lower raw float bits changed");
		check(Float.floatToRawIntBits(FamilyGateCodec.readUpper(lp,1))==bits,"Upper raw float bits changed");
		check(Float.floatToRawIntBits(FamilyGateCodec.readFinalBlosum(lp,1))==bits,"Final BLOSUM raw float bits changed");
	}

	private static void directed(double value){
		finalPair(value,value);
		final float low=FamilyGateCodec.lower(value);
		check(Float.isFinite(low) && (double)low<=value,"Downward lower inequality");
		check((double)Math.nextUp(low)>value,"Lower must be greatest representable float not above original");
		if(value>-Float.MAX_VALUE){
			final float high=FamilyGateCodec.roundUpper(value);
			check(Float.isFinite(high) && (double)high>=value,"Upward upper inequality");
			check((double)Math.nextDown(high)<value,"Upper must be least representable float not below original");
		}
	}

	private static void finalPair(double lower,double upper){
		final double bound=Math.max(lower,Math.max(0,Math.max(upper*1.1,upper+0.005)));
		if(bound>Float.MAX_VALUE){reject(()->FamilyGateCodec.finalBlosum(lower,upper)); return;}
		final float value=FamilyGateCodec.finalBlosum(lower,upper);
		check(Float.isFinite(value) && (double)value>=bound,"Final bound must cover the raw uncapped-margin formula");
		check((double)Math.nextDown(value)<bound,"Final bound must be least float not below the raw formula");
		final byte[] encoded=FamilyGateCodec.appendFinalBlosum(new ByteBuilder(),value).toBytes();
		final double decoded=FamilyGateCodec.readFinalBlosum(parser(encoded),0);
		check(decoded==value,"Final cutoff float bits must survive roundtrip");
		check(!(bound>decoded) && !(Math.nextDown(decoded)>decoded),"Raw formula and below-cutoff score must reject");
		check(!(value>value) && Math.nextUp(value)>value,"Strict float comparison rejects equality and accepts next float");
	}

	private static LineParser1 parser(byte[] token){return new LineParser1('\t').set(token);}
	private static void reject(Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){checks++; return;}
		throw new AssertionError("Malformed codec input was accepted");
	}
	private static void check(boolean condition,String message){if(!condition){throw new AssertionError(message);} checks++;}
	private static int checks;
}
