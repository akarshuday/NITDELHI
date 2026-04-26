  <!-- Features Carousel Script -->
  <script>

    /* â”€â”€ Features Carousel Logic â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
    /* â”€â”€ Features Carousel Logic (Self-Contained) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
    (function () {
      let currentSlide = 0;
      const slideCount = 4;
      let carouselInterval;
      let typeInterval;

      window.initCarousel = function () {
        console.log("Carousel: Initializing...");
        const dots = document.querySelectorAll('.nav-dot');
        const track = document.getElementById('carousel-track');
        const titles = document.querySelectorAll('.feature-title');

        // Store original text for typewriter effect
        titles.forEach(t => t.dataset.text = t.textContent);

        if (!track) return;

        dots.forEach(dot => {
          dot.addEventListener('click', () => {
            const index = parseInt(dot.dataset.index);
            goToSlide(index);
            resetCarouselTimer();
          });
        });

        startCarouselTimer();
        // Trigger first typewriter
        if (titles[0]) typeWriter(titles[0], titles[0].dataset.text);
        console.log("Carousel: Ready.");
      };

      function typeWriter(el, text) {
        if (typeInterval) clearInterval(typeInterval);
        el.textContent = '';
        let i = 0;
        typeInterval = setInterval(() => {
          if (i < text.length) {
            el.textContent += text.charAt(i);
            i++;
          } else {
            clearInterval(typeInterval);
          }
        }, 75); // Slower typing speed (75ms)
      }

      function goToSlide(index) {
        currentSlide = index;
        const track = document.getElementById('carousel-track');
        const dots = document.querySelectorAll('.nav-dot');
        if (!track) return;

        track.style.transform = `translateX(-${currentSlide * 100}%)`;
        dots.forEach((dot, i) => {
          dot.classList.toggle('active', i === currentSlide);
        });

        // Re-trigger typewriter for current slide
        const activeSlide = track.children[index];
        if (activeSlide) {
          const title = activeSlide.querySelector('.feature-title');
          if (title) typeWriter(title, title.dataset.text);
        }
      }

      function nextSlide() {
        currentSlide = (currentSlide + 1) % slideCount;
        goToSlide(currentSlide);
      }

      function startCarouselTimer() {
        if (carouselInterval) clearInterval(carouselInterval);
        carouselInterval = setInterval(nextSlide, 5000); // Slower interval (5s)
      }

      function resetCarouselTimer() {
        startCarouselTimer();
      }
    })();
