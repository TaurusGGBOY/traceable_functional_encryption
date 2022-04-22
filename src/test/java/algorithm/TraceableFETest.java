package algorithm;


import it.unisa.dia.gas.jpbc.Element;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.*;

public class TraceableFETest {
    @Test
    public void testTraceableFE() throws IOException {
        int lambda = 256;
        int k = 10;
        int n = 5;
        int id = 2;

        TraceableFE traceableFe = new TraceableFE();

        // setup
        LocalDateTime startSetup = LocalDateTime.now();
        TraceableFE.SetupParams setupParams = traceableFe.setup(lambda, k, n);
        LocalDateTime endSetup = LocalDateTime.now();
        TraceableFE.Params params = setupParams.PK.params;

        // init x y
        List<Element> x = new ArrayList<>();
        List<Element> y = new ArrayList<>();
        Element xySum = params.pairing.getZr().newZeroElement();
        for (int i = 0; i < k; i++) {
            x.add(params.pairing.getZr().newRandomElement().getImmutable());
            y.add(params.pairing.getZr().newRandomElement().getImmutable());
            xySum = xySum.add(x.get(i).mul(y.get(i)));
        }

        // extract
        LocalDateTime startKeyDer = LocalDateTime.now();
        TraceableFE.SKX SKX = traceableFe.extract(params, id, setupParams.msk, x);
        LocalDateTime endKeyDer = LocalDateTime.now();

        // encrypt
        LocalDateTime startEncrypt = LocalDateTime.now();
        TraceableFE.CT ct = traceableFe.encrypt(setupParams.PK, y);
        LocalDateTime endEncrypt = LocalDateTime.now();

        // decrypt
        LocalDateTime startDecrypt = LocalDateTime.now();
        TraceableFE.DecryptResult decryptResult = traceableFe.decrypt(params, ct, SKX, x, id);
        LocalDateTime endDecrypt = LocalDateTime.now();

        // assert if correct
        Element eggxy = params.egg.powZn(xySum);
        assert decryptResult.decryptResult.equals(eggxy);

        // time
        Duration setupDuration = Duration.between(startSetup, endSetup);
        Duration extractDuration = Duration.between(startKeyDer, endKeyDer);
        Duration encryptDuration = Duration.between(startEncrypt, endEncrypt);
        Duration decryptDuration = Duration.between(startDecrypt, endDecrypt);

        System.out.println("setupTime: " + setupDuration.toMillis());
        System.out.println("extractTime: " + extractDuration.toMillis());
        System.out.println("encryptTime: " + encryptDuration.toMillis());
        System.out.println("decryptTime: " + decryptDuration.toMillis());
    }

}