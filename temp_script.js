    function showModal(id) {
      console.log("UI: showModal", id);
      closeModals();
      const modal = document.getElementById('modal-' + id);
      if (modal) modal.classList.add('active');
    }

    function closeModals() {
      document.querySelectorAll('.modal-overlay').forEach(m => m.classList.remove('active'));
      document.querySelectorAll('.auth-error').forEach(e => e.style.display = 'none');
    }

    function checkAuthStatus() {
      console.log("Auth: Checking status...");
      const token = localStorage.getItem('auth_token');
      const name = localStorage.getItem('user_name');
      const guestNav = document.getElementById('guest-nav');
      const userNav = document.getElementById('user-nav');
      const userDisplay = document.getElementById('user-display-name');

      if (token && name) {
        document.body.classList.remove('locked');
        if (guestNav) guestNav.style.display = 'none';
        if (userNav) userNav.style.display = 'flex';
        if (userDisplay) userDisplay.textContent = name;
        const joinContainer = document.getElementById('hero-join-container');
        if (joinContainer) joinContainer.style.display = 'none';
      } else {
        document.body.classList.add('locked');
        if (guestNav) guestNav.style.display = 'flex';
        if (userNav) userNav.style.display = 'none';
        const joinContainer = document.getElementById('hero-join-container');
        if (joinContainer) joinContainer.style.display = 'block';
      }
    }

    function logout() {
      localStorage.removeItem('auth_token');
      localStorage.removeItem('user_name');
      checkAuthStatus();
      window.location.reload();
    }

    function openLab() {
      const token = localStorage.getItem('auth_token');
      if (token) {
        showHistory();
      } else {
        showModal('login');
      }
    }

    async function showHistory() {
      showModal('history');
      const container = document.getElementById('history-container');
      if (container) container.innerHTML = '<div style="padding:40px; text-align:center; color:var(--muted)">Loading your history...</div>';
      
      try {
        const res = await fetch('/history');
        const data = await res.json();
        
        if (data.error) throw new Error(data.error);
        
        if (!data.length) {
          if (container) container.innerHTML = '<div style="padding:40px; text-align:center; color:var(--muted)">No past detections found in the lab records.</div>';
          return;
        }
        
        if (container) container.innerHTML = data.map(item => `
          <div class="history-item">
            <div class="history-meta">${new Date(item.timestamp).toLocaleString()}</div>
            <div class="history-title">Target: ${item.bacteria || 'Custom Genome'}</div>
            <div class="history-results">
              Detected Markers: ${item.found_genes || 'None'}
            </div>
          </div>
        `).join('');
      } catch (e) {
        if (container) container.innerHTML = `<div style="padding:40px; text-align:center; color:#f87171">Failed to load history: ${e.message}</div>`;
      }
    }

    async function handleAuthSubmit(e, type) {
      e.preventDefault();
      const errorDiv = document.getElementById(type + '-error') || document.getElementById('forgot-error');
      if (errorDiv) errorDiv.style.display = 'none';

      let url = '/' + type;
      let body = {};

      if (type === 'register') {
        body = {
          name: document.getElementById('signup-name').value,
          email: document.getElementById('signup-email').value,
          password: document.getElementById('signup-password').value
        };
      } else if (type === 'login') {
        body = {
          email: document.getElementById('login-email').value,
          password: document.getElementById('login-password').value
        };
      } else if (type === 'forgot') {
        url = '/forgot-password';
        body = { email: document.getElementById('forgot-email').value };
      } else if (type === 'reset') {
        url = '/reset-password';
        body = {
          email: document.getElementById('forgot-email').value,
          token: document.getElementById('reset-otp').value,
          new_password: document.getElementById('reset-password').value
        };
      }

      try {
        const resp = await fetch(url, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(body)
        });
        const data = await resp.json();
        if (!resp.ok) {
          if (errorDiv) {
            errorDiv.textContent = data.error || 'Request failed';
            errorDiv.style.display = 'block';
          }
          return;
        }

        if (type === 'register') {
          alert('Registration successful! Please log in.');
          showModal('login');
        } else if (type === 'login') {
          localStorage.setItem('auth_token', data.token);
          localStorage.setItem('user_name', data.name);
          closeModals();
          checkAuthStatus();
          window.location.reload();
        } else if (type === 'forgot') {
          document.getElementById('forgot-form').style.display = 'none';
          document.getElementById('reset-form').style.display = 'block';
          document.getElementById('forgot-success').style.display = 'block';
          if (data.otp) {
            alert("Your OTP is: " + data.otp);
          }
        } else if (type === 'reset') {
          alert('Password reset successful! Please log in.');
          showModal('login');
          document.getElementById('forgot-form').style.display = 'block';
          document.getElementById('reset-form').style.display = 'none';
          document.getElementById('forgot-success').style.display = 'none';
        }
      } catch (err) {
        if (errorDiv) {
          errorDiv.textContent = 'Connection error';
          errorDiv.style.display = 'block';
        }
      }
    }

    // Initial check
