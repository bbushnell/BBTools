package server;

import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;

import shared.Parse;
import shared.Timer;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;

public class StressTest2 {
    public static void main(String[] args) {
        long iterations = (args.length > 0 ? Parse.parseKMG(args[0]) : 10000);
        String payload = (args.length > 1 ? args[1] : "NC_012345");
        String message = "https://taxonomy.jgi.doe.gov/accession/" + payload;
        
        int maxConcurrent = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(maxConcurrent);
        
        final HttpClient client = HttpClient.newBuilder()
            .version(HttpClient.Version.HTTP_2)
            .build();
        
        List<Future<Void>> futures = new ArrayList<Future<Void>>();
        
        Timer t = new Timer();
        for (long i = 0; i < iterations; i++) {
            final HttpRequest request = HttpRequest.newBuilder().uri(URI.create(message)).build();
            futures.add(executor.submit(new Callable<Void>() {
                @Override
                public Void call() throws Exception {
                    try {
                        client.send(request, HttpResponse.BodyHandlers.ofString());
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    return null;
                }
            }));
        }
        
        // Wait for all futures to complete
        for (Future<Void> future : futures) {
            try {
                future.get();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        
        executor.shutdown();
        t.stop("Sent " + iterations + " queries in ");
    }
}
