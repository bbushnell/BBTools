package test.server;

import java.io.IOException;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.time.Duration;

import parse.Parse;
import shared.Timer;

/**
 * Serial Java 11 HttpClient load helper for finite local endpoints.
 * Uses 10 second connect and 30 second request timeouts so a silent fixture
 * endpoint cannot hang the test indefinitely.
 */
public class StressTest3 {

	private static final Duration CONNECT_TIMEOUT=Duration.ofSeconds(10);
	private static final Duration REQUEST_TIMEOUT=Duration.ofSeconds(30);

	/**
	 * Sends blocking GET requests serially to an explicit endpoint.
	 * @param args endpoint, optional iteration count
	 */
	public static void main(String[] args){
		if(args.length<1 || args.length>2){
			throw new IllegalArgumentException("Usage: java test.server.StressTest3 <endpoint> [iterations=10000]");
		}
		URI endpoint=URI.create(args[0]);
		long iterations=(args.length>1 ? Parse.parseKMG(args[1]) : 10000);
		if(iterations<1){throw new IllegalArgumentException("Iterations must be positive: "+iterations);}

		HttpClient client=HttpClient.newBuilder()
				.version(HttpClient.Version.HTTP_2)
				.connectTimeout(CONNECT_TIMEOUT)
				.build();
		HttpRequest request=HttpRequest.newBuilder()
				.uri(endpoint)
				.timeout(REQUEST_TIMEOUT)
				.build();
		Stats stats=new Stats();
		Timer t=new Timer();

		for(long i=0; i<iterations; i++){
			try{
				stats.record(client.send(request, HttpResponse.BodyHandlers.discarding()).statusCode());
			}catch(InterruptedException e){
				Thread.currentThread().interrupt();
				stats.exceptions++;
				break;
			}catch(IOException e){
				stats.exceptions++;
			}
		}

		t.stop("Completed "+stats.completed+" of "+iterations+" serial queries in ");
		stats.finish(iterations);
	}

	/** Response counters and failure policy independent of Java assertions. */
	private static class Stats {
		void record(int status){
			completed++;
			if(status<200 || status>299){non2xx++;}
		}
		void finish(long expected){
			System.out.println("completed\t"+completed);
			System.out.println("non2xx\t"+non2xx);
			System.out.println("exceptions\t"+exceptions);
			if(completed!=expected || non2xx>0 || exceptions>0){
				throw new IllegalStateException("completed="+completed+", expected="+expected+
						", non2xx="+non2xx+", exceptions="+exceptions);
			}
		}
		long completed;
		long non2xx;
		long exceptions;
	}
}
